"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib
import sys

import flow
import numpy as np
from flow import environments


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


# ____________________________________________________________________________
"""Setting progress label"""


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_created_box(job):
    """Check if the lammps simulation box has been created for the job."""
    return job.isfile("box.lammps") and job.isfile("box.json")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_copy_files(job):
    """Check if the submission scripts have been copied over for the job."""
    return job.isfile("submit.pbs") and job.isfile("in.minimize")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_minimized_equilibrated_nvt(job):
    """Check if the lammps minimization step has run for the job."""
    return job.isfile("minimized.restart-0")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@flow.with_job
def lammps_equilibrated_npt(job):
    """Check if the lammps equilibration step has run and passed is_equilibrated for the job."""
    import pathlib

    import numpy as np

    from reproducibility_project.src.analysis.equilibration import (
        is_equilibrated,
    )

    p = pathlib.Path(".")
    list_of_filenames = list(p.glob("eqlog*.txt"))
    # grab the filename with the largest number
    counter = -1
    latest_eqdata = False
    for file in list_of_filenames:
        step = int(file.name[5:].split(".")[0])
        if step > counter:
            counter = step
            latest_eqdata = file
    if latest_eqdata:
        data = np.genfromtxt(latest_eqdata.name, skip_header=1)
        check_equil = [
            is_equilibrated(data[:, 1])[0],
            is_equilibrated(data[:, 2])[0],
            is_equilibrated(data[:, 4])[0],
            is_equilibrated(data[:, 6])[0],
        ]
    else:
        check_equil = False
    return job.isfile("equilibrated-npt.restart") and np.all(check_equil)


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_production_npt(job):
    """Check if the lammps production step has run for the job."""
    return job.isfile("production-npt.restart")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_production_nvt(job):
    """Check if the lammps nvt production step has run for the job."""
    return job.isfile("production-nvt.restart")


# sample job to get decorrelated data


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_reformatted_data(job):
    """Check if lammps has output density information for the job."""
    return job.isfile("log-npt.txt") and job.isfile("log-nvt.txt")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_created_gsd(job):
    """Check if the mdtraj has converted the production to a gsd trajectory for the job."""
    return job.isfile("trajectory-npt.gsd") and job.isfile("trajectory-nvt.gsd")


# _____________________________________________________________________
"""Setting up workflow operation"""


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.post(lammps_created_box)
@flow.with_job
def built_lammps(job):
    """Create initial configurations of the system statepoint."""
    import foyer
    from mbuild.formats.lammpsdata import write_lammpsdata

    from reproducibility_project.src.molecules.system_builder import (
        construct_system,
    )
    from reproducibility_project.src.utils.forcefields import load_ff

    system = construct_system(job.sp)[0]
    parmed_structure = system.to_parmed()
    ff = load_ff(job.sp.forcefield_name)
    system.save(
        "box.json"
    )  # save the compound as a json object for reading back in to mbuild
    typed_box = ff.apply(parmed_structure)
    typed_box.save(
        "box.top"
    )  # save to gromacs topology for later conversions in mdtraj
    typed_box.save(
        "box.gro"
    )  # save to gromacs topology for later conversions in mdtraj
    write_lammpsdata(
        typed_box,
        "box.lammps",
        atom_style="full",
        unit_style="real",
        mins=[system.get_boundingbox().vectors[0]],
        maxs=[system.get_boundingbox().vectors[1]],
        use_rb_torsions=True,
    )  # write out a lammps topology
    return


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(lammps_created_box)
@Project.post(lammps_copy_files)
@flow.with_job
@flow.cmd
def lammps_cp_files(job):
    """Copy over run files for lammps and the PBS scheduler."""
    lmps_submit_path = "../../src/engine_input/lammps-VU/submit.pbs"
    lmps_run_path = "../../src/engine_input/lammps-VU/in.*"
    msg = f"cp {lmps_submit_path} {lmps_run_path} ./"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(lammps_copy_files)
@Project.post(lammps_minimized_equilibrated_nvt)
@flow.with_job
@flow.cmd
def lammps_em_nvt(job):
    """Run energy minimization and nvt ensemble."""
    if job.sp.molecule == "ethanolAA":
        tstep = 1.0
    else:
        tstep = 2.0
    in_script_name = "in.minimize"
    r_cut = job.sp.r_cut * 10
    if job.sp.long_range_correction:
        pass_lrc = "yes"
    else:
        pass_lrc = "no"
    if job.sp.cutoff_style == "shift":
        pass_shift = "yes"
    else:
        pass_shift = "no"
    modify_submit_scripts(in_script_name, job.id)
    msg = f"qsub -v 'infile={in_script_name}, seed={job.sp.replica+1}, T={job.sp.temperature}, P={job.sp.pressure}, rcut={r_cut}, tstep={tstep} lrc={pass_lrc} shift={pass_shift}' submit.pbs"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(lammps_minimized_equilibrated_nvt)
@Project.post(lammps_equilibrated_npt)
@flow.with_job
@flow.cmd
def lammps_equil_npt(job):
    """Run npt ensemble equilibration."""
    if job.sp.molecule == "ethanolAA":
        tstep = 1.0
    else:
        tstep = 2.0
    in_script_name = "in.equilibration"
    modify_submit_scripts(in_script_name, job.id)
    r_cut = job.sp.r_cut * 10
    msg = f"qsub -v 'infile={in_script_name}, seed={job.sp.replica+1}, T={job.sp.temperature}, P={job.sp.pressure}, rcut={r_cut}, tstep={tstep}' submit.pbs"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(lammps_equilibrated_npt)
@Project.post(lammps_production_npt)
@flow.with_job
@flow.cmd
def lammps_prod_npt(job):
    """Run npt ensemble production."""
    if job.sp.molecule == "ethanolAA":
        tstep = 1.0
    else:
        tstep = 2.0
    in_script_name = "in.production-npt"
    modify_submit_scripts(in_script_name, job.id)
    r_cut = job.sp.r_cut * 10
    msg = f"qsub -v 'infile={in_script_name}, seed={job.sp.replica+1}, T={job.sp.temperature}, P={job.sp.pressure}, rcut={r_cut}, tstep={tstep}' submit.pbs"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(lammps_production_npt)
@Project.post(lammps_production_nvt)
@flow.with_job
@flow.cmd
def lammps_prod_nvt(job):
    """Run npt ensemble production."""
    if job.sp.molecule == "ethanolAA":
        tstep = 1.0
    else:
        tstep = 2.0
    in_script_name = "in.production-nvt"
    modify_submit_scripts(in_script_name, job.id)
    r_cut = job.sp.r_cut * 10
    msg = f"qsub -v 'infile={in_script_name}, seed={job.sp.replica+1}, T={job.sp.temperature}, P={job.sp.pressure}, rcut={r_cut}, tstep={tstep}' submit.pbs"

    return msg


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(lammps_production_nvt)
@Project.post(lammps_reformatted_data)
@flow.with_job
def lammps_reformat_data(job):
    """Take data from thermo.txt and reformat to log.txt with correct units.

    Lammps units real: energy=kcal/mol, temp=K, press=atm, density=g/cm^3, step=2fs
    Project units: energy=kJ/mol, temp=K, press=MPa, density=amu/nm^3, step=1ps
    """
    import numpy as np
    import pandas as pd

    df_in = pd.read_csv(job.ws + "/prlog-npt.txt", delimiter=" ", header=0)
    attr_list = ["step", "pe", "ke", "press", "temp", "density"]
    new_titles_list = [
        "timestep",
        "potential_energy",
        "kinetic_energy",
        "pressure",
        "temperature",
        "density",
    ]
    # convert units
    KCAL_TO_KJ = 4.184  # kcal to kj
    ATM_TO_MPA = 0.101325  # atm to mpa
    GPCM3_TO_AMUPNM3 = 0.6023  # g/cm^3 to amu/nm^3
    df_in["pe"] = df_in["pe"] * KCAL_TO_KJ
    df_in["ke"] = df_in["ke"] * KCAL_TO_KJ
    df_in["press"] = df_in["press"] * ATM_TO_MPA
    df_in["density"] = df_in["density"] * GPCM3_TO_AMUPNM3
    df_out = df_in[attr_list]
    df_out.columns = new_titles_list
    df_out.to_csv("log-npt.txt", header=True, index=False, sep=" ")

    df_in = pd.read_csv(job.ws + "/prlog-nvt.txt", delimiter=" ", header=0)
    attr_list = ["step", "pe", "ke", "press", "temp", "density"]
    new_titles_list = [
        "timestep",
        "potential_energy",
        "kinetic_energy",
        "pressure",
        "temperature",
        "density",
    ]
    df_in["pe"] = df_in["pe"] * KCAL_TO_KJ
    df_in["ke"] = df_in["ke"] * KCAL_TO_KJ
    df_in["press"] = df_in["press"] * ATM_TO_MPA
    df_in["density"] = df_in["density"] * GPCM3_TO_AMUPNM3
    df_out = df_in[attr_list]
    df_out.columns = new_titles_list
    df_out.to_csv("log-nvt.txt", header=True, index=False, sep=" ")


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(lammps_reformatted_data)
@Project.post(lammps_created_gsd)
@flow.with_job
def lammps_create_gsd(job):
    """Create an rdf from the gsd file using Freud analysis scripts."""
    # Create rdf data from the production run
    import mdtraj as md

    traj = md.load("prod-npt.xtc", top="box.gro")
    traj.save("trajectory-npt.gsd")
    traj = md.load("prod-nvt.xtc", top="box.gro")
    traj.save("trajectory-nvt.gsd")
    return


def modify_submit_scripts(filename, jobid, cores=8):
    """Modify the submission scripts to include the job and simulation type in the header."""
    with open("submit.pbs", "r") as f:
        lines = f.readlines()
        lines[1] = "#PBS -N {}-{}\n".format(filename[3:], jobid[0:4])
    with open("submit.pbs", "w") as f:
        f.writelines(lines)
    return


if __name__ == "__main__":
    pr = Project()
    pr.main()
