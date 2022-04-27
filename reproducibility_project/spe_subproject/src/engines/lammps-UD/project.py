"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib

import flow
import foyer
import numpy as np
from flow import environments
from mbuild.formats.lammpsdata import write_lammpsdata

from reproducibility_project.src.utils.forcefields import load_ff


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


# ____________________________________________________________________________
"""Setting progress label"""


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-UD")
def CreatedEngineInput(job):
    """Check if the .json molecule topology was converted to engine input."""
    return job.isfile("in.spe") and job.isfile("box.lammps")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-UD")
def OutputThermoData(job):
    """Check if the engine loaded the input files and wrote out thermo data."""
    return job.isfile("prlog-npt.txt")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-UD")
def FinishedSPECalc(job):
    """Check if the log-spe.txt has been created."""
    return job.isfile("log-spe.txt")


# _____________________________________________________________________
"""Setting up workflow operation"""


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-UD")
@Project.post(CreatedEngineInput)
@flow.with_job
@flow.cmd
def LoadSystemSnapShot(job):
    """Create initial configurations of the system statepoint."""
    import mbuild as mb

    pr = Project()
    snapshot_directory = (
        pathlib.Path(pr.root_directory()) / "src" / "system_snapshots"
    )
    molecule = job.sp.molecule
    molecule_filename = molecule + ".json"
    box = mb.load(str(snapshot_directory / molecule_filename))
    parmed_box = box.to_parmed()
    ff = load_ff(job.sp.forcefield_name)
    # Apply forcefield and write out engine input files
    # __________________________________________________
    typed_box = ff.apply(parmed_box)
    write_lammpsdata(
        typed_box,
        "box.lammps",
        atom_style="full",
        unit_style="real",
        mins=[box.get_boundingbox().vectors[0]],
        maxs=[box.get_boundingbox().vectors[1]],
        use_rb_torsions=True,
    )  # write out a lammps topology
    lmps_submit_path = "../../src/engine_input/lammps-UD/submit.slurm"
    lmps_run_path = "../../src/engine_input/lammps-UD/in.*"
    msg = f"cp {lmps_submit_path} {lmps_run_path} ./"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-UD")
@Project.pre(CreatedEngineInput)
@Project.post(OutputThermoData)
@flow.with_job
@flow.cmd
def CalculateEnergy(job):
    """Load onto a cluster and output the point energy for the snapshot."""
    # __________________________________________________
    tstep = 2.0
    in_script_name = "in.spe"
    modify_submit_scripts(in_script_name, job.id)
    if job.sp.molecule in ["waterSPCE"]:
        add_shake(in_script_name, 14)
    if job.sp.molecule in ["waterSPCE", "ethanolAA"]:
        add_pppm(in_script_name, 12)
        modify_engine_scripts(
            in_script_name, 7, "pair_style lj/cut/coul/long ${rcut}\n"
        )
    if job.sp.molecule in ["ethanolAA"]:
        add_14coul(in_script_name, 28)
    r_cut = job.sp.r_cut * 10
    pass_lrc = "yes"
    pass_shift = "no"
    msg = f"sbatch submit.slurm {in_script_name} {job.sp.replica+1} {job.sp.temperature} {job.sp.pressure} {r_cut} {tstep} {pass_lrc} {pass_shift}"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-UD")
@Project.pre(OutputThermoData)
@Project.post(FinishedSPECalc)
@flow.with_job
@flow.cmd
def FormatTextFile(job):
    """Take the output from the simulation engine and convert it to log-spe.txt for data comparisons.

    See README.md for spe_subproject for formatting information.
    """
    # __________________________________________________
    import numpy as np
    import pandas as pd

    attr_list = [
        "step",
        "total",
        "pot",
        "vdw",
        "coul",
        "pair",
        "bonds",
        "angles",
        "dihedrals",
        "tail",
        "kspace",
    ]
    df_npt_in = pd.read_csv(
        job.ws + "/prlog-npt.txt", delimiter=" ", comment="#", names=attr_list
    )
    print(df_npt_in)
    KCAL_TO_KJ = 4.184  # kcal to kj
    ATM_TO_MPA = 0.101325  # atm to mpa
    for attr in attr_list:
        df_npt_in[attr] = df_npt_in[attr] * KCAL_TO_KJ
    df_npt_in.to_csv("log-spe.txt", header=True, index=False, sep=" ")


def add_shake(filename, ln):
    with open(filename, "r") as f:
        lines = f.readlines()
        lines[ln] = "fix fix_shake all shake 0.00001 20 1000 b 1 a 1\n"
    with open(filename, "w") as f:
        f.writelines(lines)
    return


def add_14coul(filename, ln):
    with open(filename, "r") as f:
        lines = f.readlines()
        lines[ln] = "special_bonds lj/coul 0 0 0.5\n"
    with open(filename, "w") as f:
        f.writelines(lines)
    return


def add_pppm(filename, ln):
    with open(filename, "r") as f:
        lines = f.readlines()
        lines[ln] = "kspace_style pppm 0.00001\n"
    with open(filename, "w") as f:
        f.writelines(lines)
    return


def remove_shake(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        lines[27] = "\n"
    with open(filename, "w") as f:
        f.writelines(lines)
    return


def modify_submit_scripts(filename, jobid, cores=8):
    """Modify the submission scripts to include the job and simulation type in the header."""
    with open("submit.slurm", "r") as f:
        lines = f.readlines()
        lines[1] = "#SBATCH --job-name={}-{}\n".format(filename[3:], jobid[0:4])
    with open("submit.slurm", "w") as f:
        f.writelines(lines)
    return


def modify_engine_scripts(filename, ln, info):
    """Modify any line of any scripts."""
    with open(filename, "r") as f:
        lines = f.readlines()
        lines[ln] = info
    with open(filename, "w") as f:
        f.writelines(lines)
    return


if __name__ == "__main__":
    pr = Project()
    pr.main()
