"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib
import sys

import flow
import numpy as np
from flow import environments
from flow.environment import DefaultSlurmEnvironment


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


class Rahman(DefaultSlurmEnvironment):
    """Subclass of DefaultPBSEnvironment for VU's Rahman cluster."""

    # template = "rahman_lmp.sh"


# ____________________________________________________________________________
"""Setting progress label"""


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def CreatedEngineInput(job):
    """Check if the .json molecule topology was converted to engine input."""
    return job.isfile("box.lammps")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def OutputThermoData(job):
    """Check if the engine loaded the input files and wrote out thermo data."""
    return job.isfile("prlog-npt.txt")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def FinishedSPECalc(job):
    """Check if the log-spe.txt has been created."""
    return job.isfile("log-spe.txt")


# _____________________________________________________________________
"""Setting up workflow operation"""


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.post(CreatedEngineInput)
@flow.with_job
def LoadSystemSnapShot(job):
    """Create initial configurations of the system statepoint."""
    import mbuild as mb
    from mbuild.formats.lammpsdata import write_lammpsdata

    from reproducibility_project.src.utils.forcefields import load_ff

    pr = Project()
    snapshot_directory = (
        pathlib.Path(pr.root_directory()) / "src" / "system_snapshots"
    )
    molecule = job.sp.molecule
    molecule_filename = molecule + ".json"
    box = mb.load(str(snapshot_directory / molecule_filename))
    parmed_structure = box.to_parmed()
    ff = load_ff(job.sp.forcefield_name)
    typed_box = ff.apply(parmed_structure)
    typed_box.save(
        "box.top", overwrite=True
    )  # save to gromacs topology for later conversions in mdtraj
    typed_box.save(
        "box.gro", overwrite=True
    )  # save to gromacs topology for later conversions in mdtraj
    write_lammpsdata(
        typed_box,
        "box.lammps",
        atom_style="full",
        unit_style="real",
        mins=[box.box.vectors[0]],
        maxs=[box.box.vectors[1]],
        use_rb_torsions=True,
    )  # write out a lammps topology


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(CreatedEngineInput)
@Project.post(OutputThermoData)
@flow.with_job
@flow.cmd
def CalculateEnergy(job):
    """Load onto a cluster and output the point energy for the snapshot.

    Copy over run files for lammps and the PBS scheduler.
    """
    lmps_submit_path = "../../src/engine_input/lammps-VU/submit.pbs"
    lmps_run_path = "../../src/engine_input/lammps-VU/in.production-npt"
    msg = f"cp {lmps_submit_path} {lmps_run_path} ./"
    os.system(msg)
    """Run energy minimization and nvt ensemble."""
    in_script_name = "submit.pbs"
    modify_submit_scripts(in_script_name, job.id)
    in_script_name = "in.production-npt"
    r_cut = job.sp.r_cut * 10
    if job.sp.molecule == "ethanolAA":
        tstep = 1.0
    else:
        tstep = 2.0

    if (
        "waterSPCE" in job.sp.molecule or "ethanolAA" in job.sp.molecule
    ):  # add charges for water and ethanol
        modify_engine_scripts(
            in_script_name, "pair_style lj/cut/coul/long ${rcut}\n", 7
        )
        modify_engine_scripts(
            in_script_name,
            "kspace_style pppm 1.0e-5 #PPPM Ewald, relative error in forces\n",
            12,
        )
    if "waterSPCE" == job.sp.molecule:  # Fix SHAKE for spce water
        modify_engine_scripts(
            in_script_name, "fix rigbod all shake 0.00001 20 0 b 1 a 1\n", 14
        )
    modify_engine_scripts(
        in_script_name, " special_bonds lj/coul 0 0 0.5\n", 16
    )  # use 1-4 combining rules for lammps
    msg = f"qsub -v 'infile={in_script_name}, seed={job.sp.replica+1}, T={job.sp.temperature}, P={job.sp.pressure}, rcut={r_cut}, tstep={tstep}' submit.pbs"

    return msg


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(OutputThermoData)
@Project.post(FinishedSPECalc)
@flow.with_job
def FormatTextFile(job):
    """Take data from thermo.txt and reformat to log.txt with correct units.

    Lammps units real: energy=kcal/mol, temp=K, press=atm, density=g/cm^3, step=2fs
    Project units: energy=kJ/mol, temp=K, press=MPa, density=g/cm^3, step=1ps
    """
    import numpy as np
    import pandas as pd

    df_in = pd.read_csv(job.ws + "/prlog-npt.txt", delimiter=" ", header=0)
    attr_list = [
        "pe",
        "evdwl",
        "ecoul",
        "epair",
        "ebond",
        "eangle",
        "edihed",
        "etail",
        "elong",
    ]
    new_titles_list = [
        "potential_energy",
        "tot_vdw_energy",
        "short_range_electrostatics",
        "pair_energy",
        "bonds_energy",
        "angles_energy",
        "dihedrals_energy",
        "tail_energy",
        "long_range_electrostatics",
    ]
    # convert units
    KCAL_TO_KJ = 4.184  # kcal to kj
    df_in = df_in * KCAL_TO_KJ
    df_out = df_in[attr_list]
    df_out.columns = new_titles_list
    # calculate new values
    df_out["tot_electrostatics"] = (
        df_out["short_range_electrostatics"]
        + df_out["long_range_electrostatics"]
    )
    df_out["tot_bonded_energy"] = (
        df_out["bonds_energy"]
        + df_out["angles_energy"]
        + df_out["dihedrals_energy"]
    )
    df_out["tot_pair_energy"] = (
        df_out["tot_vdw_energy"] + df_out["tot_electrostatics"]
    )
    df_out["intramolecular_energy"] = None
    df_out["intermolecular_energy"] = None
    df_out.to_csv("log-spe.txt", header=True, index=False, sep=",")


def modify_submit_scripts(filename, jobid, cores=8):
    """Modify the submission scripts to include the job and simulation type in the header."""
    with open("submit.pbs", "r") as f:
        lines = f.readlines()
        lines[1] = "#PBS -N {}\n".format(jobid)
    with open("submit.pbs", "w") as f:
        f.writelines(lines)


def modify_engine_scripts(filename, msg, line):
    """Modify the submission scripts to include the job and simulation type in the header."""
    with open(filename, "r") as f:
        lines = f.readlines()
        lines[line] = msg
    with open(filename, "w") as f:
        f.writelines(lines)


if __name__ == "__main__":
    pr = Project()
    pr.main()
