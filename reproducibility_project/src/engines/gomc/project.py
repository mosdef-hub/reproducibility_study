"""GOMC's setup for signac, signac-flow, signac-dashboard for this study."""
# project.py

import os
import subprocess

import flow
import matplotlib.pyplot as plt

# from flow.environment import StandardEnvironment
import mbuild as mb
import mbuild.formats.charmm_writer as mf_charmm
import mbuild.formats.gomc_conf_writer as gomc_control
import numpy as np
import pandas as pd
import pymbar
import signac
import unyt as u
from flow import FlowProject
from flow.environment import DefaultSlurmEnvironment

from reproducibility_project.src.analysis.equilibration import is_equilibrated
from reproducibility_project.src.molecules.system_builder import (
    construct_system,
)
from reproducibility_project.src.utils.forcefields import get_ff_path
from reproducibility_project.src.utils.plotting import plot_data_with_t0_line


class Project(FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


class Grid(DefaultSlurmEnvironment):  # Grid(StandardEnvironment):
    """Subclass of DefaultSlurmEnvironment for WSU's Grid cluster."""

    hostname_pattern = r".*\.grid\.wayne\.edu"
    template = "grid.sh"


# ******************************************************
# users typical variables, but not all (start)
# ******************************************************
# set binary path to gomc binary files (the bin folder).
# If the gomc binary files are callable directly from the terminal without a path,
# please just enter and empty string (i.e., "" or '')
gomc_binary_path = "/Users/brad/Programs/GOMC/GOMC_dev_8_25_21/bin"

# number of MC cycles
MC_cycles_melt_equilb_NVT = 5 * 10 ** 3  # set value for paper = 5 * 10 ** 3
MC_cycles_equilb_NVT = 5 * 10 ** 3  # set value for paper = 5 * 10 ** 3
MC_cycles_equilb_design_ensemble = (
    40 * 10 ** 3
)  # set value for paper = 40 * 10 ** 3
MC_cycles_production = 120 * 10 ** 3  # set value for paper = 120 * 10 ** 3

output_data_every_X_MC_cycles = 10

# max number of equilibrium selected runs
equilb_design_ensemble_max_number = 3

# force field (FF) file for all simulations in that job
# Note: do not add extensions
ff_filename_str = "in_FF"

# initial mosdef structure and coordinates
# Note: do not add extensions
mosdef_structure_box_0_name_str = "mosdef_box_0"
mosdef_structure_box_1_name_str = "mosdef_box_1"

# melt equilb simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
melt_equilb_NVT_control_file_name_str = "melt_NVT"

# equilb simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
equilb_NVT_control_file_name_str = "equilb_NVT"

# The equilb using the ensemble used for the simulation design, which
# includes the simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
equilb_design_ensemble_control_file_name_str = "equilb_design_ensemble"

# The production run using the ensemble used for the simulation design, which
# includes the simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
production_control_file_name_str = "production_run"


path_from_job_to_box_inputs = "../../"

walltime_mosdef_hr = 4
walltime_gomc_hr = 368
memory_needed = 16

use_pymbar = True  # True of False

ff_info_dict = {
    "trappe-ua": {
        "ngpu": 0,
        "ncpu": 2,
        "Ewald": False,
        "ElectroStatic": False,
        "VDWGeometricSigma": False,
    },
    "benzene-ua": {
        "ngpu": 0,
        "ncpu": 2,
        "Ewald": False,
        "ElectroStatic": False,
        "VDWGeometricSigma": False,
    },
    "spce": {
        "ngpu": 1,
        "ncpu": 4,
        "Ewald": True,
        "ElectroStatic": True,
        "VDWGeometricSigma": False,
    },
    "oplsaa": {
        "ngpu": 1,
        "ncpu": 4,
        "Ewald": True,
        "ElectroStatic": True,
        "VDWGeometricSigma": True,
    },
}

# ******************************************************
# users typical variables, but not all (end)
# ******************************************************


# ******************************************************
# signac and GOMC-MOSDEF code (start)
# ******************************************************

# define project
project = signac.init_project("mosdef_reproducibility")

# ******************************************************
# ******************************************************
# create some initial variable to be store in each jobs
# directory in an additional json file, and test
# to see if they are written (start).
# ******************************************************
# ******************************************************
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
def part_1a_initial_data_input_to_json(job):
    """Check that the initial job data is written to the json files."""
    data_written_bool = False
    if job.isfile(f"{'signac_job_document.json'}"):
        data_written_bool = True

    return data_written_bool


@Project.pre(lambda j: j.sp.engine == "gomc")
@Project.post(part_1a_initial_data_input_to_json)
@Project.operation.with_directives(
    {
        "np": lambda job: ff_info_dict.get(job.sp.forcefield_name).get("ncpu"),
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@flow.with_job
def initial_parameters(job):
    """Set the initial job parameters into the jobs doc json file."""
    # select
    job.doc.ngpu = ff_info_dict.get(job.sp.forcefield_name).get("ngpu")
    if job.doc.ngpu == 0:
        job.doc.cpu_or_gpu = "CPU"
    elif job.doc.ngpu == 1:
        job.doc.cpu_or_gpu = "GPU"
    else:
        raise ValueError(
            "CPU and GPU can not be deterimined as force field (FF) is not available in the selection."
        )

    # FF type to directory and path
    job.doc.forcefield_directory_name = get_ff_path(job.sp.forcefield_name)

    # reformat ensembles for input to GOMC
    if job.sp.ensemble in ["GEMC-NPT", "GEMC_NPT"]:
        job.doc.production_ensemble = "GEMC_NPT"
    elif job.sp.ensemble in ["GEMC-NVT", "GEMC_NVT"]:
        job.doc.production_ensemble = "GEMC_NVT"
    else:
        job.doc.production_ensemble = job.sp.ensemble

    # list replica seed numbers
    replica_no_to_seed_dict = {
        0: 0,
        1: 1,
        2: 2,
        3: 3,
        4: 4,
        5: 5,
        6: 6,
        7: 7,
        8: 8,
        9: 9,
        10: 10,
        11: 11,
        12: 12,
        13: 13,
        14: 14,
        15: 15,
        16: 16,
        17: 17,
        18: 18,
        19: 19,
        20: 20,
    }

    job.doc.replica_number_int = replica_no_to_seed_dict.get(
        int(job.sp.replica)
    )

    # set rcut, ewalds
    job.doc.Ewald = ff_info_dict.get(job.sp.forcefield_name).get("Ewald")
    job.doc.ElectroStatic = ff_info_dict.get(job.sp.forcefield_name).get(
        "ElectroStatic"
    )
    job.doc.VDWGeometricSigma = ff_info_dict.get(job.sp.forcefield_name).get(
        "VDWGeometricSigma"
    )

    # set the initial iteration number of the simulation
    job.doc.equilb_design_ensemble_dict = {}
    job.doc.equilb_design_ensemble_number = 0
    job.doc.equilb_design_ensemble_max_number = (
        equilb_design_ensemble_max_number
    )
    job.doc.equilb_design_ensemble_max_number_under_limit = True
    job.doc.stable_equilb_design_ensemble = False

    job.doc.production_run_ensemble_dict = {}

    if job.doc.production_ensemble in ["NVT", "NPT"]:
        job.doc.N_liquid = job.sp.N_liquid
        job.doc.N_vap = 0
        job.doc.N_total = job.sp.N_liquid
    elif job.doc.production_ensemble in ["GCMC", "GEMC_NVT", "GEMC_NPT"]:
        job.doc.N_liquid = job.sp.N_liquid
        job.doc.N_vap = job.sp.N_vap
        job.doc.N_total = job.sp.N_liquid + job.sp.N_vap
    # reject if set design ensemble is NVT
    # currently pymbar is done of of density, which will not work for NVT
    if job.doc.production_ensemble == "NVT":
        raise ValueError(
            "ERROR: The NVT ensemble is not currently available for this project.py "
            "script, as pymbar is done based on density, which will not work for NVT"
        )

    # select binary path and binary file
    job.doc.gomc_binary_path = gomc_binary_path

    if job.doc.production_ensemble in ["NPT", "NVT"]:
        job.doc.melt_NVT_gomc_binary_file = f"GOMC_{job.doc.cpu_or_gpu}_NVT"
        job.doc.equilb_NVT_gomc_binary_file = f"GOMC_{job.doc.cpu_or_gpu}_NVT"
        job.doc.equilb_design_ensemble_gomc_binary_file = (
            f"GOMC_{job.doc.cpu_or_gpu}_NPT"
        )
    elif job.doc.production_ensemble in ["GEMC_NVT", "GEMC_NPT"]:
        job.doc.melt_NVT_gomc_binary_file = f"GOMC_{job.doc.cpu_or_gpu}_GEMC"
        job.doc.equilb_NVT_gomc_binary_file = f"GOMC_{job.doc.cpu_or_gpu}_GEMC"
        job.doc.equilb_design_ensemble_gomc_binary_file = (
            f"GOMC_{job.doc.cpu_or_gpu}_GEMC"
        )
    elif job.doc.production_ensemble in ["GCMC"]:
        job.doc.melt_NVT_gomc_binary_file = f"GOMC_{job.doc.cpu_or_gpu}_GCMC"
        job.doc.equilb_NVT_gomc_binary_file = f"GOMC_{job.doc.cpu_or_gpu}_GCMC"
        job.doc.equilb_design_ensemble_gomc_binary_file = (
            f"GOMC_{job.doc.cpu_or_gpu}_GCMC"
        )
    else:
        raise ValueError(
            "ERROR: A wrong ensemble has been specified for the gomc binary file"
        )

    if job.doc.production_ensemble in ["NPT"]:
        job.doc.production_ensemble_gomc_binary_file = (
            f"GOMC_{job.doc.cpu_or_gpu}_NPT"
        )
    elif job.doc.production_ensemble in ["NVT"]:
        job.doc.production_ensemble_gomc_binary_file = (
            f"GOMC_{job.doc.cpu_or_gpu}_NVT"
        )
    elif job.doc.production_ensemble in ["GEMC_NVT", "GEMC_NPT"]:
        job.doc.production_ensemble_gomc_binary_file = (
            f"GOMC_{job.doc.cpu_or_gpu}_GEMC"
        )
    elif job.doc.production_ensemble in ["GCMC"]:
        job.doc.production_ensemble_gomc_binary_file = (
            f"GOMC_{job.doc.cpu_or_gpu}_GCMC"
        )
    else:
        raise ValueError(
            "ERROR: A wrong ensemble has been specified for the gomc binary file"
        )


# ******************************************************
# ******************************************************
# create some initial variable to be store in each jobs
# directory in an additional json file, and test
# to see if they are written (end).
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC psf, pdb, and force field (FF) files were written (start)
# ******************************************************
# ******************************************************


@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
def part_1b_under_equilb_design_ensemble_run_limit(job):
    """Check that the equilbrium design ensemble run is under it's run limit."""
    try:
        if (
            job.doc.equilb_design_ensemble_number
            >= job.doc.equilb_design_ensemble_max_number
        ):
            job.doc.equilb_design_ensemble_max_number_under_limit = False
            return job.doc.equilb_design_ensemble_max_number_under_limit

        else:
            return True
    except:
        return False


# check if GOMC-MOSDEF wrote the gomc files
# @Project.pre(select_production_ensemble)
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def mosdef_input_written(job):
    """Check that the mosdef files (psf, pdb, and force field (FF) files) are written ."""
    file_written_bool = False

    if job.doc.production_ensemble in ["NPT", "NVT"]:
        if (
            job.isfile(f"{path_from_job_to_box_inputs}/{ff_filename_str}.inp")
            and job.isfile(
                f"{path_from_job_to_box_inputs}/{mosdef_structure_box_0_name_str}.psf"
            )
            and job.isfile(
                f"{path_from_job_to_box_inputs}/{mosdef_structure_box_0_name_str}.pdb"
            )
        ):
            file_written_bool = True
    elif job.doc.production_ensemble in ["GCMC", "GEMC_NPT", "GEMC_NPT"]:
        if (
            job.isfile(f"{path_from_job_to_box_inputs}/{ff_filename_str}.inp")
            and job.isfile(
                f"{path_from_job_to_box_inputs}/{mosdef_structure_box_0_name_str}.psf"
            )
            and job.isfile(
                f"{path_from_job_to_box_inputs}/{mosdef_structure_box_0_name_str}.pdb"
            )
            and job.isfile(
                f"{path_from_job_to_box_inputs}/{mosdef_structure_box_1_name_str}.psf"
            )
            and job.isfile(
                f"{path_from_job_to_box_inputs}/{mosdef_structure_box_1_name_str}.pdb"
            )
        ):
            file_written_bool = True
    return file_written_bool


# ******************************************************
# ******************************************************
# check if GOMC psf, pdb, and FF files were written (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC control file was written (start)
# ******************************************************
# ******************************************************
# function for checking if the GOMC control file is written
def gomc_control_file_written(job, control_filename_str):
    """General check that the gomc control files are written."""
    file_written_bool = False
    control_file = f"{control_filename_str}.conf"

    if job.isfile(control_file):
        with open(f"{control_file}", "r") as fp:
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "OutputName" in line:
                    split_move_line = line.split()
                    if split_move_line[0] == "OutputName":
                        file_written_bool = True

    return file_written_bool


# checking if the GOMC control file is written for the melt equilb NVT run
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_2a_melt_equilb_NVT_control_file_written(job):
    """General check that the melt_equilb_NVT_control (high temperature) gomc control file is written."""
    return gomc_control_file_written(job, melt_equilb_NVT_control_file_name_str)


# checking if the GOMC control file is written for the equilb NVT run
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_2b_equilb_NVT_control_file_written(job):
    """General check that the equilb_NVT_control (run temperature) gomc control file is written."""
    return gomc_control_file_written(job, equilb_NVT_control_file_name_str)


# checking if the GOMC control file is written for the equilb run with the selected ensemble
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_2c_equilb_design_ensemble_control_file_written(job):
    """General check that the equilb_design_ensemble (run temperature) gomc control file is written."""
    try:
        if job.doc.equilb_design_ensemble_max_number_under_limit is True:
            return gomc_control_file_written(
                job,
                job.doc.equilb_design_ensemble_dict[
                    str(job.doc.equilb_design_ensemble_number)
                ]["output_name_control_file_name"],
            )
    except:
        return False


# checking if the GOMC control file is written for the production run
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_2d_production_control_file_written(job):
    """General check that the production run (run temperature) gomc control file is written."""
    try:
        if job.doc.equilb_design_ensemble_max_number_under_limit is True:
            return gomc_control_file_written(
                job,
                job.doc.production_run_ensemble_dict[
                    str(job.doc.equilb_design_ensemble_number)
                ]["input_name_control_file_name"],
            )
    except:
        return False


# ******************************************************
# ******************************************************
# check if GOMC control file was written (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC simulations started (start)
# ******************************************************
# ******************************************************
# function for checking if GOMC simulations are started
def gomc_simulation_started(job, control_filename_str):
    """General check to see if the gomc simulation is started."""
    output_started_bool = False
    if job.isfile("out_{}.dat".format(control_filename_str)) and job.isfile(
        "{}_merged.psf".format(control_filename_str)
    ):
        output_started_bool = True

    return output_started_bool


# check if melt equilb_NVT GOMC run is started by seeing if the GOMC consol file and the merged psf exist
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_3a_output_melt_equilb_NVT_started(job):
    """Check to see if the melt_equilb_NVT (high temperature) gomc simulation is started."""
    return gomc_simulation_started(job, melt_equilb_NVT_control_file_name_str)


# check if equilb_NVT GOMC run is started by seeing if the GOMC consol file and the merged psf exist
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_3b_output_equilb_NVT_started(job):
    """Check to see if the equilb_NVT (set temperature) gomc simulation is started."""
    return gomc_simulation_started(job, equilb_NVT_control_file_name_str)


# check if equilb_with design ensemble GOMC run is started by seeing if the GOMC consol file and the merged psf exist
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_3c_output_equilb_design_ensemble_started(job):
    """Check to see if the equilb_design_ensemble (set temperature) gomc simulation is started."""
    if job.isfile(
        "out_{}.dat".format(
            job.doc.equilb_design_ensemble_dict[
                str(job.doc.equilb_design_ensemble_number)
            ]["output_name_control_file_name"]
        )
    ):
        if job.doc.equilb_design_ensemble_max_number_under_limit is True:
            return gomc_simulation_started(
                job,
                job.doc.equilb_design_ensemble_dict[
                    str(job.doc.equilb_design_ensemble_number)
                ]["output_name_control_file_name"],
            )

    else:
        return False


# check if production GOMC run is started by seeing if the GOMC consol file and the merged psf exist
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_3d_output_production_run_started(job):
    """Check to see if the production run (set temperature) gomc simulation is started."""
    if job.isfile(
        "out_{}.dat".format(
            job.doc.production_run_ensemble_dict[
                str(job.doc.equilb_design_ensemble_number)
            ]["output_name_control_file_name"]
        )
    ):
        if job.doc.equilb_design_ensemble_max_number_under_limit is True:
            return gomc_simulation_started(
                job,
                job.doc.production_run_ensemble_dict[
                    str(job.doc.equilb_design_ensemble_number)
                ]["output_name_control_file_name"],
            )
    else:
        return False


# ******************************************************
# ******************************************************
# check if GOMC simulations started (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC simulation are completed properly (start)
# ******************************************************
# ******************************************************

# function for checking if GOMC simulations are completed properly
def gomc_sim_completed_properly(job, control_filename_str):
    """General check to see if the gomc simulation was completed properly."""
    job_run_properly_bool = False
    output_log_file = "out_{}.dat".format(control_filename_str)
    if job.isfile(output_log_file):
        # with open(f"workspace/{job.id}/{output_log_file}", "r") as fp:
        with open(f"{output_log_file}", "r") as fp:
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "Move" in line:
                    split_move_line = line.split()
                    if (
                        split_move_line[0] == "Move"
                        and split_move_line[1] == "Type"
                        and split_move_line[2] == "Mol."
                        and split_move_line[3] == "Kind"
                    ):
                        job_run_properly_bool = True
    else:
        job_run_properly_bool = False

    return job_run_properly_bool


# check if melt equilb NVT GOMC run completed by checking the end of the GOMC consol file
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_4a_job_melt_equilb_NVT_completed_properly(job):
    """Check to see if the melt_equilb_NVT (high temperature) gomc simulation was completed properly."""
    return gomc_sim_completed_properly(
        job, melt_equilb_NVT_control_file_name_str
    )


# check if equilb NVT GOMC run completed by checking the end of the GOMC consol file
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_4b_job_equilb_NVT_completed_properly(job):
    """Check to see if the equilb_NVT (set temperature) gomc simulation was completed properly."""
    return gomc_sim_completed_properly(job, equilb_NVT_control_file_name_str)


# check if equilb selected ensemble GOMC run completed by checking the end of the GOMC consol file
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_4c_job_equilb_design_ensemble_completed_properly(job):
    """Check to see if the equilb_design_ensemble (set temperature) gomc simulation was completed properly."""
    try:
        return gomc_sim_completed_properly(
            job,
            job.doc.equilb_design_ensemble_dict[
                str(job.doc.equilb_design_ensemble_number)
            ]["output_name_control_file_name"],
        )

    except:
        return False


# check if production GOMC run completed by checking the end of the GOMC consol file
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_4d_job_production_run_completed_properly(job):
    """Check to see if the production run (set temperature) gomc simulation was completed properly."""
    return gomc_sim_completed_properly(job, production_control_file_name_str)


# ******************************************************
# ******************************************************
# check if GOMC simulation are completed properly (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# build system, with option to write the force field (force field (FF)), pdb, psf files.
# Note: this is needed to write GOMC control file, even if a restart (start)
# ******************************************************
# build system
def build_charmm(job, write_files=True):
    """Build the Charmm object and potentially write the pdb, psd, and force field (FF) files."""
    print("#**********************")
    print("Started: GOMC Charmm Object")
    print("#**********************")

    print("Running: box packing")
    [box_0, box_1] = construct_system(job.sp)
    print("Completed: box packing")

    Bead_to_atom_name_dict = {
        "_CH4": "C",
        "_CH3": "C",
        "_CH2": "C",
        "_CH": "C",
        "_HC": "C",
    }
    Molecule_ResName_List = [job.sp.molecule]

    if job.sp.molecule in ["waterSPCE"]:
        gomc_fix_bonds_angles_list = Molecule_ResName_List
    else:
        gomc_fix_bonds_angles_list = None

    if job.doc.production_ensemble in ["NVT", "NPT"]:
        charmm = mf_charmm.Charmm(
            box_0,
            mosdef_structure_box_0_name_str,
            structure_box_1=None,
            filename_box_1=None,
            ff_filename=ff_filename_str,
            forcefield_selection=job.doc.forcefield_directory_name,
            residues=Molecule_ResName_List,
            bead_to_atom_name_dict=Bead_to_atom_name_dict,
            gomc_fix_bonds_angles=gomc_fix_bonds_angles_list,
        )

    elif job.doc.production_ensemble in ["GCMC", "GEMC_NVT", "GEMC_NPT"]:
        charmm = mf_charmm.Charmm(
            box_0,
            mosdef_structure_box_0_name_str,
            structure_box_1=box_1,
            filename_box_1=mosdef_structure_box_1_name_str,
            ff_filename=ff_filename_str,
            forcefield_selection=job.doc.forcefield_directory_name,
            residues=Molecule_ResName_List,
            bead_to_atom_name_dict=Bead_to_atom_name_dict,
            gomc_fix_bonds_angles=gomc_fix_bonds_angles_list,
        )

    if write_files == True:
        charmm.write_inp()

        charmm.write_psf()

        charmm.write_pdb()

    print("#**********************")
    print("Completed: GOMC Charmm Object")
    print("#**********************")

    return charmm


# ******************************************************
# ******************************************************
# build system, with option to write the force field (FF), pdb, psf files.
# Note: this is needed to write GOMC control file, even if a restart (end)
# ******************************************************

# ******************************************************
# ******************************************************
# Creating GOMC files (pdb, psf, force field (FF), and gomc control files (start)
# ******************************************************
# ******************************************************
@Project.pre(lambda j: j.sp.engine == "gomc")
@Project.pre(part_1a_initial_data_input_to_json)
@Project.pre(part_1b_under_equilb_design_ensemble_run_limit)
@Project.post(part_2a_melt_equilb_NVT_control_file_written)
@Project.post(part_2b_equilb_NVT_control_file_written)
@Project.post(part_2c_equilb_design_ensemble_control_file_written)
@Project.post(part_2d_production_control_file_written)
@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@flow.with_job
def build_psf_pdb_ff_gomc_conf(job):
    """Build the Charmm object and write the pdb, psd, and force field (FF) files for all the simulations in the workspace."""
    charmm_object_with_files = build_charmm(job, write_files=True)

    # ******************************************************
    # melt_NVT - psf, pdb, force field (FF) file writing and GOMC control file writing  (start)
    # ******************************************************
    print("#**********************")
    print("Starting: melt_NVT GOMC control file writing")
    print("#**********************")

    Restart = False

    temperature = (1000 * u.K).to_value("K")
    pressure = (job.sp.pressure * u.kPa).to_value("bar")

    output_name_control_file_name = melt_equilb_NVT_control_file_name_str
    restart_control_file_name_str = None

    # replica unmber cycles
    MC_cycles = MC_cycles_melt_equilb_NVT
    # calc MC steps

    MC_steps = int(MC_cycles * job.doc.N_total)
    EqSteps = 1000

    seed_no = job.doc.replica_number_int

    production_ensemble = job.doc.production_ensemble

    Ewald = ff_info_dict.get(job.sp.forcefield_name).get("Ewald")
    ElectroStatic = ff_info_dict.get(job.sp.forcefield_name).get(
        "ElectroStatic"
    )
    VDWGeometricSigma = ff_info_dict.get(job.sp.forcefield_name).get(
        "VDWGeometricSigma"
    )
    Rcut = job.sp.r_cut * u.nm
    Rcut = Rcut.to_value("angstrom")
    if job.sp.cutoff_style == "hard":
        LRC = False
    else:
        raise ValueError("ERROR: Not a valid cutoff_style")

    # output all data and calc frequecy
    output_true_list_input = [
        True,
        int(output_data_every_X_MC_cycles * job.doc.N_total),
    ]
    output_false_list_input = [
        False,
        int(output_data_every_X_MC_cycles * job.doc.N_total),
    ]

    if production_ensemble in ["NVT", "NPT"]:
        used_ensemble = "NVT"

        if job.sp.molecule in ["methaneUA"]:
            VolFreq = (None,)
            SwapFreq = (None,)
            DisFreq = (1.00,)
            RotFreq = (None,)
            RegrowthFreq = (None,)

        elif job.sp.molecule in ["pentaneUA", "benzeneUA", "ethanolAA"]:
            VolFreq = (None,)
            SwapFreq = (None,)
            DisFreq = (0.34,)
            RotFreq = (0.33,)
            RegrowthFreq = (0.33,)

        elif job.sp.molecule in ["waterSPCE"]:
            VolFreq = (None,)
            SwapFreq = (None,)
            DisFreq = (0.5,)
            RotFreq = (0.5,)
            RegrowthFreq = (None,)

        else:
            raise ValueError(
                "Moleules MC move rations not listed in the GOMC control file writer."
            )

        if Restart is True:
            Coordinates_box_0 = "{}_BOX_0_restart.pdb".format(
                output_name_control_file_name
            )
            Structure_box_0 = "{}_BOX_0_restart.psf".format(
                output_name_control_file_name
            )
            binCoordinates_box_0 = "{}_BOX_0_restart.coor".format(
                output_name_control_file_name
            )
            extendedSystem_box_0 = "{}_BOX_0_restart.xsc".format(
                output_name_control_file_name
            )

        elif Restart is False:
            Coordinates_box_0 = None
            Structure_box_0 = None
            binCoordinates_box_0 = None
            extendedSystem_box_0 = None

        Coordinates_box_1 = None
        Structure_box_1 = None
        binCoordinates_box_1 = None
        extendedSystem_box_1 = None

    elif job.doc.production_ensemble in ["GCMC", "GEMC_NVT", "GEMC_NPT"]:
        used_ensemble = job.doc.production_ensemble

        if job.sp.molecule in ["methaneUA"]:
            VolFreq = (None,)
            SwapFreq = (None,)
            DisFreq = (1.00,)
            RotFreq = (None,)
            RegrowthFreq = (None,)

        elif job.sp.molecule in ["pentaneUA", "benzeneUA", "ethanolAA"]:
            VolFreq = (None,)
            SwapFreq = (None,)
            DisFreq = (0.34,)
            RotFreq = (0.33,)
            RegrowthFreq = (0.33,)

        elif job.sp.molecule in ["waterSPCE"]:
            VolFreq = (None,)
            SwapFreq = (None,)
            DisFreq = (0.5,)
            RotFreq = (0.5,)
            RegrowthFreq = (None,)

        else:
            raise ValueError(
                "Moleules MC move rations not listed in the GOMC control file writer."
            )

        if Restart is True:
            Coordinates_box_0 = "{}_BOX_0_restart.pdb".format(
                output_name_control_file_name
            )
            Structure_box_0 = "{}_BOX_0_restart.psf".format(
                output_name_control_file_name
            )
            binCoordinates_box_0 = "{}_BOX_0_restart.coor".format(
                output_name_control_file_name
            )
            extendedSystem_box_0 = "{}_BOX_0_restart.xsc".format(
                output_name_control_file_name
            )

            Coordinates_box_1 = "{}_BOX_1_restart.pdb".format(
                output_name_control_file_name
            )
            Structure_box_1 = "{}_BOX_1_restart.psf".format(
                output_name_control_file_name
            )
            binCoordinates_box_1 = "{}_BOX_1_restart.coor".format(
                output_name_control_file_name
            )
            extendedSystem_box_1 = "{}_BOX_1_restart.xsc".format(
                output_name_control_file_name
            )

        elif Restart is False:
            Coordinates_box_0 = None
            Structure_box_0 = None
            binCoordinates_box_0 = None
            extendedSystem_box_0 = None

            Coordinates_box_1 = None
            Structure_box_1 = None
            binCoordinates_box_1 = None
            extendedSystem_box_1 = None

    gomc_control.write_gomc_control_file(
        charmm_object_with_files,
        output_name_control_file_name,
        used_ensemble,
        MC_steps,
        temperature,
        ff_psf_pdb_file_directory=None,
        check_input_files_exist=False,
        Parameters="{}.inp".format(ff_filename_str),
        Restart=Restart,
        RestartCheckpoint=False,
        ExpertMode=True,
        Coordinates_box_0=Coordinates_box_0,
        Structure_box_0=Structure_box_0,
        binCoordinates_box_0=binCoordinates_box_0,
        extendedSystem_box_0=extendedSystem_box_0,
        binVelocities_box_0=None,
        Coordinates_box_1=Coordinates_box_1,
        Structure_box_1=Structure_box_1,
        binCoordinates_box_1=binCoordinates_box_1,
        extendedSystem_box_1=extendedSystem_box_1,
        binVelocities_box_1=None,
        input_variables_dict={
            "PRNG": seed_no,
            "Pressure": pressure,
            "Ewald": Ewald,
            "ElectroStatic": ElectroStatic,
            "VDWGeometricSigma": VDWGeometricSigma,
            "Rcut": Rcut,
            "VolFreq": VolFreq[-1],
            "SwapFreq": SwapFreq[-1],
            "DisFreq": DisFreq[-1],
            "RotFreq": RotFreq[-1],
            "RegrowthFreq": RegrowthFreq[-1],
            "OutputName": output_name_control_file_name,
            "EqSteps": EqSteps,
            "PressureCalc": output_true_list_input,
            "RestartFreq": output_true_list_input,
            "CheckpointFreq": output_true_list_input,
            "ConsoleFreq": output_true_list_input,
            "BlockAverageFreq": output_true_list_input,
            "HistogramFreq": output_true_list_input,
            "CoordinatesFreq": output_false_list_input,
            "DCDFreq": output_false_list_input,
            "LRC": LRC,
            "RcutLow": 1,
            "CBMC_First": 12,
            "CBMC_Nth": 10,
            "CBMC_Ang": 50,
            "CBMC_Dih": 50,
        },
    )

    print("#**********************")
    print("Completed: melt_NVT GOMC control file writing")
    print("#**********************")
    # ******************************************************
    # melt_NVT - psf, pdb, force field (FF) file writing and GOMC control file writing  (end)
    # ******************************************************

    # ******************************************************
    # equilb_NVT - GOMC control file writing  (start)
    # ******************************************************
    print("#**********************")
    print("Started: equilb_NVT GOMC control file writing")
    print("#**********************")

    Restart = True

    temperature = (job.sp.temperature * u.K).to_value("K")
    pressure = (job.sp.pressure * u.kPa).to_value("bar")

    output_name_control_file_name = equilb_NVT_control_file_name_str
    restart_control_file_name_str = melt_equilb_NVT_control_file_name_str

    # replica unmber cycles
    MC_cycles = MC_cycles_equilb_NVT
    # calc MC steps
    MC_steps = int(MC_cycles * job.doc.N_total)
    EqSteps = 1000

    seed_no = job.doc.replica_number_int

    production_ensemble = job.doc.production_ensemble

    Ewald = ff_info_dict.get(job.sp.forcefield_name).get("Ewald")
    ElectroStatic = ff_info_dict.get(job.sp.forcefield_name).get(
        "ElectroStatic"
    )
    VDWGeometricSigma = ff_info_dict.get(job.sp.forcefield_name).get(
        "VDWGeometricSigma"
    )
    Rcut = job.sp.r_cut * u.nm
    Rcut = Rcut.to_value("angstrom")
    if job.sp.cutoff_style == "hard":
        LRC = False
    else:
        raise ValueError("ERROR: Not a valid cutoff_style")

    # output all data and calc frequecy
    output_true_list_input = [
        True,
        int(output_data_every_X_MC_cycles * job.doc.N_total),
    ]
    output_false_list_input = [
        False,
        int(output_data_every_X_MC_cycles * job.doc.N_total),
    ]

    if production_ensemble in ["NVT", "NPT"]:
        used_ensemble = "NVT"

        if job.sp.molecule in ["methaneUA"]:
            VolFreq = (None,)
            SwapFreq = (None,)
            DisFreq = (1.00,)
            RotFreq = (None,)
            RegrowthFreq = (None,)

        elif job.sp.molecule in ["pentaneUA", "benzeneUA", "ethanolAA"]:
            VolFreq = (None,)
            SwapFreq = (None,)
            DisFreq = (0.34,)
            RotFreq = (0.33,)
            RegrowthFreq = (0.33,)

        elif job.sp.molecule in ["waterSPCE"]:
            VolFreq = (None,)
            SwapFreq = (None,)
            DisFreq = (0.5,)
            RotFreq = (0.5,)
            RegrowthFreq = (None,)

        else:
            raise ValueError(
                "Moleules MC move rations not listed in the GOMC control file writer."
            )

        if Restart is True:
            Coordinates_box_0 = "{}_BOX_0_restart.pdb".format(
                restart_control_file_name_str
            )
            Structure_box_0 = "{}_BOX_0_restart.psf".format(
                restart_control_file_name_str
            )
            binCoordinates_box_0 = "{}_BOX_0_restart.coor".format(
                restart_control_file_name_str
            )
            extendedSystem_box_0 = "{}_BOX_0_restart.xsc".format(
                restart_control_file_name_str
            )

        elif Restart is False:
            Coordinates_box_0 = None
            Structure_box_0 = None
            binCoordinates_box_0 = None
            extendedSystem_box_0 = None

        Coordinates_box_1 = None
        Structure_box_1 = None
        binCoordinates_box_1 = None
        extendedSystem_box_1 = None

    elif job.doc.production_ensemble in ["GCMC", "GEMC_NVT", "GEMC_NPT"]:
        used_ensemble = job.doc.production_ensemble

        if job.sp.molecule in ["methaneUA"]:
            VolFreq = (None,)
            SwapFreq = (None,)
            DisFreq = (1.00,)
            RotFreq = (None,)
            RegrowthFreq = (None,)

        elif job.sp.molecule in ["pentaneUA", "benzeneUA", "ethanolAA"]:
            VolFreq = (None,)
            SwapFreq = (None,)
            DisFreq = (0.34,)
            RotFreq = (0.33,)
            RegrowthFreq = (0.33,)

        elif job.sp.molecule in ["waterSPCE"]:
            VolFreq = (None,)
            SwapFreq = (None,)
            DisFreq = (0.5,)
            RotFreq = (0.5,)
            RegrowthFreq = (None,)

        else:
            raise ValueError(
                "Moleules MC move rations not listed in the GOMC control file writer."
            )

        if Restart is True:
            Coordinates_box_0 = "{}_BOX_0_restart.pdb".format(
                restart_control_file_name_str
            )
            Structure_box_0 = "{}_BOX_0_restart.psf".format(
                restart_control_file_name_str
            )
            binCoordinates_box_0 = "{}_BOX_0_restart.coor".format(
                restart_control_file_name_str
            )
            extendedSystem_box_0 = "{}_BOX_0_restart.xsc".format(
                restart_control_file_name_str
            )

            Coordinates_box_1 = "{}_BOX_1_restart.pdb".format(
                restart_control_file_name_str
            )
            Structure_box_1 = "{}_BOX_1_restart.psf".format(
                restart_control_file_name_str
            )
            binCoordinates_box_1 = "{}_BOX_1_restart.coor".format(
                restart_control_file_name_str
            )
            extendedSystem_box_1 = "{}_BOX_1_restart.xsc".format(
                restart_control_file_name_str
            )

        elif Restart is False:
            Coordinates_box_0 = None
            Structure_box_0 = None
            binCoordinates_box_0 = None
            extendedSystem_box_0 = None

            Coordinates_box_1 = None
            Structure_box_1 = None
            binCoordinates_box_1 = None
            extendedSystem_box_1 = None

    gomc_control.write_gomc_control_file(
        charmm_object_with_files,
        output_name_control_file_name,
        used_ensemble,
        MC_steps,
        temperature,
        ff_psf_pdb_file_directory=None,
        check_input_files_exist=False,
        Parameters="{}.inp".format(ff_filename_str),
        Restart=Restart,
        RestartCheckpoint=False,
        ExpertMode=True,
        Coordinates_box_0=Coordinates_box_0,
        Structure_box_0=Structure_box_0,
        binCoordinates_box_0=binCoordinates_box_0,
        extendedSystem_box_0=extendedSystem_box_0,
        binVelocities_box_0=None,
        Coordinates_box_1=Coordinates_box_1,
        Structure_box_1=Structure_box_1,
        binCoordinates_box_1=binCoordinates_box_1,
        extendedSystem_box_1=extendedSystem_box_1,
        binVelocities_box_1=None,
        input_variables_dict={
            "PRNG": seed_no,
            "Pressure": pressure,
            "Ewald": Ewald,
            "ElectroStatic": ElectroStatic,
            "VDWGeometricSigma": VDWGeometricSigma,
            "Rcut": Rcut,
            "VolFreq": VolFreq[-1],
            "SwapFreq": SwapFreq[-1],
            "DisFreq": DisFreq[-1],
            "RotFreq": RotFreq[-1],
            "RegrowthFreq": RegrowthFreq[-1],
            "OutputName": output_name_control_file_name,
            "EqSteps": EqSteps,
            "PressureCalc": output_true_list_input,
            "RestartFreq": output_true_list_input,
            "CheckpointFreq": output_true_list_input,
            "ConsoleFreq": output_true_list_input,
            "BlockAverageFreq": output_true_list_input,
            "HistogramFreq": output_true_list_input,
            "CoordinatesFreq": output_false_list_input,
            "DCDFreq": output_false_list_input,
            "LRC": LRC,
            "RcutLow": 1,
            "CBMC_First": 12,
            "CBMC_Nth": 10,
            "CBMC_Ang": 50,
            "CBMC_Dih": 50,
        },
    )

    print("#**********************")
    print("Completed: equilb_NVT GOMC control file writing")
    print("#**********************")
    # ******************************************************
    # equilb_NVT - GOMC control file writing  (end)
    # ******************************************************

    # ******************************************************
    # equilb selected_ensemble, if NVT -> NPT - GOMC control file writing  (start)
    # Note: the control files are written for the max number of equilb_design_ensemble runs
    # so the Charmm object only needs created 1 time.
    # ******************************************************
    print("#**********************")
    print("Started: equilb NPT or GEMC-NVT GOMC control file writing")
    print("#**********************")

    for number_sims_i in range(0, job.doc.equilb_design_ensemble_max_number):
        Restart = True

        temperature = (job.sp.temperature * u.K).to_value("K")
        pressure = (job.sp.pressure * u.kPa).to_value("bar")

        if number_sims_i == 0:
            restart_control_file_name_str = equilb_NVT_control_file_name_str
            output_name_control_file_name = "{}_number_{}".format(
                equilb_design_ensemble_control_file_name_str, number_sims_i
            )

        elif number_sims_i <= job.doc.equilb_design_ensemble_max_number:
            restart_control_file_name_str = "{}_number_{}".format(
                equilb_design_ensemble_control_file_name_str,
                int(number_sims_i - 1),
            )

            output_name_control_file_name = "{}_number_{}".format(
                equilb_design_ensemble_control_file_name_str, number_sims_i
            )
        job.doc.equilb_design_ensemble_dict.update(
            {
                number_sims_i: {
                    "restart_control_file_name": restart_control_file_name_str,
                    "output_name_control_file_name": output_name_control_file_name,
                }
            }
        )

        # replica unmber cycles
        MC_cycles = MC_cycles_equilb_design_ensemble

        # calc MC steps
        MC_steps = int(MC_cycles * job.doc.N_total)
        EqSteps = 1000

        seed_no = job.doc.replica_number_int

        production_ensemble = job.doc.production_ensemble

        Ewald = ff_info_dict.get(job.sp.forcefield_name).get("Ewald")
        ElectroStatic = ff_info_dict.get(job.sp.forcefield_name).get(
            "ElectroStatic"
        )
        VDWGeometricSigma = ff_info_dict.get(job.sp.forcefield_name).get(
            "VDWGeometricSigma"
        )
        Rcut = job.sp.r_cut * u.nm
        Rcut = Rcut.to_value("angstrom")
        if job.sp.cutoff_style == "hard":
            LRC = False
        else:
            raise ValueError("ERROR: Not a valid cutoff_style")

        # output all data and calc frequecy
        output_true_list_input = [
            True,
            int(output_data_every_X_MC_cycles * job.doc.N_total),
        ]
        output_false_list_input = [
            False,
            int(output_data_every_X_MC_cycles * job.doc.N_total),
        ]

        if production_ensemble in ["NVT", "NPT"]:
            used_ensemble = "NPT"

            if job.sp.molecule in ["methaneUA"]:
                VolFreq = (0.01,)
                SwapFreq = (None,)
                DisFreq = (0.99,)
                RotFreq = (None,)
                RegrowthFreq = (None,)

            elif job.sp.molecule in ["pentaneUA", "benzeneUA", "ethanolAA"]:
                VolFreq = (0.01,)
                SwapFreq = (None,)
                DisFreq = (0.33,)
                RotFreq = (0.33,)
                RegrowthFreq = (0.33,)

            elif job.sp.molecule in ["waterSPCE"]:
                VolFreq = (0.01,)
                SwapFreq = (None,)
                DisFreq = (0.49,)
                RotFreq = (0.5,)
                RegrowthFreq = (None,)

            else:
                raise ValueError(
                    "Moleules MC move ratios not listed in the GOMC control file writer."
                )

            if Restart is True:
                Coordinates_box_0 = "{}_BOX_0_restart.pdb".format(
                    restart_control_file_name_str
                )
                Structure_box_0 = "{}_BOX_0_restart.psf".format(
                    restart_control_file_name_str
                )
                binCoordinates_box_0 = "{}_BOX_0_restart.coor".format(
                    restart_control_file_name_str
                )
                extendedSystem_box_0 = "{}_BOX_0_restart.xsc".format(
                    restart_control_file_name_str
                )

            elif Restart is False:
                Coordinates_box_0 = None
                Structure_box_0 = None
                binCoordinates_box_0 = None
                extendedSystem_box_0 = None

            Coordinates_box_1 = None
            Structure_box_1 = None
            binCoordinates_box_1 = None
            extendedSystem_box_1 = None

        elif job.doc.production_ensemble in ["GCMC", "GEMC_NVT", "GEMC_NPT"]:
            used_ensemble = job.doc.production_ensemble

            if job.sp.molecule in ["methaneUA"]:
                VolFreq = (0.01,)
                SwapFreq = (0.29,)
                DisFreq = (0.70,)
                RotFreq = (None,)
                RegrowthFreq = (None,)

            elif job.sp.molecule in ["pentaneUA", "benzeneUA", "ethanolAA"]:
                VolFreq = (0.01,)
                SwapFreq = (0.20,)
                DisFreq = (0.27,)
                RotFreq = (0.26,)
                RegrowthFreq = (0.26,)

            elif job.sp.molecule in ["waterSPCE"]:
                VolFreq = (0.01,)
                SwapFreq = (0.20,)
                DisFreq = (0.40,)
                RotFreq = (0.39,)
                RegrowthFreq = (None,)

            else:
                raise ValueError(
                    "Moleules MC move ratios not listed in the GOMC control file writer."
                )

            if Restart is True:
                Coordinates_box_0 = "{}_BOX_0_restart.pdb".format(
                    restart_control_file_name_str
                )
                Structure_box_0 = "{}_BOX_0_restart.psf".format(
                    restart_control_file_name_str
                )
                binCoordinates_box_0 = "{}_BOX_0_restart.coor".format(
                    restart_control_file_name_str
                )
                extendedSystem_box_0 = "{}_BOX_0_restart.xsc".format(
                    restart_control_file_name_str
                )

                Coordinates_box_1 = "{}_BOX_1_restart.pdb".format(
                    restart_control_file_name_str
                )
                Structure_box_1 = "{}_BOX_1_restart.psf".format(
                    restart_control_file_name_str
                )
                binCoordinates_box_1 = "{}_BOX_1_restart.coor".format(
                    restart_control_file_name_str
                )
                extendedSystem_box_1 = "{}_BOX_1_restart.xsc".format(
                    restart_control_file_name_str
                )

            elif Restart is False:
                Coordinates_box_0 = None
                Structure_box_0 = None
                binCoordinates_box_0 = None
                extendedSystem_box_0 = None

                Coordinates_box_1 = None
                Structure_box_1 = None
                binCoordinates_box_1 = None
                extendedSystem_box_1 = None

        gomc_control.write_gomc_control_file(
            charmm_object_with_files,
            output_name_control_file_name,
            used_ensemble,
            MC_steps,
            temperature,
            ff_psf_pdb_file_directory=None,
            check_input_files_exist=False,
            Parameters="{}.inp".format(ff_filename_str),
            Restart=Restart,
            RestartCheckpoint=False,
            ExpertMode=True,
            Coordinates_box_0=Coordinates_box_0,
            Structure_box_0=Structure_box_0,
            binCoordinates_box_0=binCoordinates_box_0,
            extendedSystem_box_0=extendedSystem_box_0,
            binVelocities_box_0=None,
            Coordinates_box_1=Coordinates_box_1,
            Structure_box_1=Structure_box_1,
            binCoordinates_box_1=binCoordinates_box_1,
            extendedSystem_box_1=extendedSystem_box_1,
            binVelocities_box_1=None,
            input_variables_dict={
                "PRNG": seed_no,
                "Pressure": pressure,
                "Ewald": Ewald,
                "ElectroStatic": ElectroStatic,
                "VDWGeometricSigma": VDWGeometricSigma,
                "Rcut": Rcut,
                "VolFreq": VolFreq[-1],
                "SwapFreq": SwapFreq[-1],
                "DisFreq": DisFreq[-1],
                "RotFreq": RotFreq[-1],
                "RegrowthFreq": RegrowthFreq[-1],
                "OutputName": output_name_control_file_name,
                "EqSteps": EqSteps,
                "PressureCalc": output_true_list_input,
                "RestartFreq": output_true_list_input,
                "CheckpointFreq": output_true_list_input,
                "ConsoleFreq": output_true_list_input,
                "BlockAverageFreq": output_true_list_input,
                "HistogramFreq": output_false_list_input,
                "CoordinatesFreq": output_false_list_input,
                "DCDFreq": output_false_list_input,
                "LRC": LRC,
                "RcutLow": 1,
                "CBMC_First": 12,
                "CBMC_Nth": 10,
                "CBMC_Ang": 50,
                "CBMC_Dih": 50,
            },
        )
        print("#**********************")
        print("Completed: equilb NPT or GEMC-NVT GOMC control file writing")
        print("#**********************")

        # ******************************************************
        # equilb selected_ensemble, if NVT -> NPT - GOMC control file writing  (end)
        # Note: the control files are written for the max number of equilb_design_ensemble runs
        # so the Charmm object only needs created 1 time.
        # ******************************************************

        # ******************************************************
        # production NPT or GEMC-NVT - GOMC control file writing  (start)
        # ******************************************************

        print("#**********************")
        print("Started: production NPT or GEMC-NVT GOMC control file writing")
        print("#**********************")

        Restart = True

        temperature = (job.sp.temperature * u.K).to_value("K")
        pressure = (job.sp.pressure * u.kPa).to_value("bar")

        input_name_control_file_name = "{}_number_{}".format(
            production_control_file_name_str, number_sims_i
        )
        output_name_control_file_name = production_control_file_name_str
        restart_control_file_name_str = "{}_number_{}".format(
            equilb_design_ensemble_control_file_name_str, int(number_sims_i)
        )
        job.doc.production_run_ensemble_dict.update(
            {
                number_sims_i: {
                    "restart_control_file_name": restart_control_file_name_str,
                    "input_name_control_file_name": input_name_control_file_name,
                    "output_name_control_file_name": output_name_control_file_name,
                }
            }
        )

        # replica unmber cycles
        MC_cycles = MC_cycles_production

        # calc MC steps
        MC_steps = int(MC_cycles * job.doc.N_total)
        EqSteps = 1000

        seed_no = job.doc.replica_number_int

        production_ensemble = job.doc.production_ensemble

        Ewald = ff_info_dict.get(job.sp.forcefield_name).get("Ewald")
        ElectroStatic = ff_info_dict.get(job.sp.forcefield_name).get(
            "ElectroStatic"
        )
        VDWGeometricSigma = ff_info_dict.get(job.sp.forcefield_name).get(
            "VDWGeometricSigma"
        )
        Rcut = job.sp.r_cut * u.nm
        Rcut = Rcut.to_value("angstrom")
        if job.sp.cutoff_style == "hard":
            LRC = False
        else:
            raise ValueError("ERROR: Not a valid cutoff_style")

        # output all data and calc frequecy
        output_true_list_input = [
            True,
            int(output_data_every_X_MC_cycles * job.doc.N_total),
        ]
        output_false_list_input = [
            False,
            int(output_data_every_X_MC_cycles * job.doc.N_total),
        ]

        if production_ensemble in ["NVT", "NPT"]:
            used_ensemble = "NPT"

            if job.sp.molecule in ["methaneUA"]:
                VolFreq = (0.01,)
                SwapFreq = (None,)
                DisFreq = (0.99,)
                RotFreq = (None,)
                RegrowthFreq = (None,)

            elif job.sp.molecule in ["pentaneUA", "benzeneUA", "ethanolAA"]:
                VolFreq = (0.01,)
                SwapFreq = (None,)
                DisFreq = (0.33,)
                RotFreq = (0.33,)
                RegrowthFreq = (0.33,)

            elif job.sp.molecule in ["waterSPCE"]:
                VolFreq = (0.01,)
                SwapFreq = (None,)
                DisFreq = (0.49,)
                RotFreq = (0.5,)
                RegrowthFreq = (None,)

            else:
                raise ValueError(
                    "Moleules MC move rations not listed in the GOMC control file writer."
                )

            if Restart is True:
                Coordinates_box_0 = "{}_BOX_0_restart.pdb".format(
                    restart_control_file_name_str
                )
                Structure_box_0 = "{}_BOX_0_restart.psf".format(
                    restart_control_file_name_str
                )
                binCoordinates_box_0 = "{}_BOX_0_restart.coor".format(
                    restart_control_file_name_str
                )
                extendedSystem_box_0 = "{}_BOX_0_restart.xsc".format(
                    restart_control_file_name_str
                )

            elif Restart is False:
                Coordinates_box_0 = None
                Structure_box_0 = None
                binCoordinates_box_0 = None
                extendedSystem_box_0 = None

            Coordinates_box_1 = None
            Structure_box_1 = None
            binCoordinates_box_1 = None
            extendedSystem_box_1 = None

        elif job.doc.production_ensemble in ["GCMC", "GEMC_NVT", "GEMC_NPT"]:
            used_ensemble = job.doc.production_ensemble

            if job.sp.molecule in ["methaneUA"]:
                VolFreq = (0.01,)
                SwapFreq = (0.29,)
                DisFreq = (0.70,)
                RotFreq = (None,)
                RegrowthFreq = (None,)

            elif job.sp.molecule in ["pentaneUA", "benzeneUA", "ethanolAA"]:
                VolFreq = (0.01,)
                SwapFreq = (0.20,)
                DisFreq = (0.27,)
                RotFreq = (0.26,)
                RegrowthFreq = (0.26,)

            elif job.sp.molecule in ["waterSPCE"]:
                VolFreq = (0.01,)
                SwapFreq = (0.20,)
                DisFreq = (0.40,)
                RotFreq = (0.39,)
                RegrowthFreq = (None,)

            else:
                raise ValueError(
                    "Moleules MC move rations not listed in the GOMC control file writer."
                )

            if Restart is True:
                Coordinates_box_0 = "{}_BOX_0_restart.pdb".format(
                    restart_control_file_name_str
                )
                Structure_box_0 = "{}_BOX_0_restart.psf".format(
                    restart_control_file_name_str
                )
                binCoordinates_box_0 = "{}_BOX_0_restart.coor".format(
                    restart_control_file_name_str
                )
                extendedSystem_box_0 = "{}_BOX_0_restart.xsc".format(
                    restart_control_file_name_str
                )

                Coordinates_box_1 = "{}_BOX_1_restart.pdb".format(
                    restart_control_file_name_str
                )
                Structure_box_1 = "{}_BOX_1_restart.psf".format(
                    restart_control_file_name_str
                )
                binCoordinates_box_1 = "{}_BOX_1_restart.coor".format(
                    restart_control_file_name_str
                )
                extendedSystem_box_1 = "{}_BOX_1_restart.xsc".format(
                    restart_control_file_name_str
                )

            elif Restart is False:
                Coordinates_box_0 = None
                Structure_box_0 = None
                binCoordinates_box_0 = None
                extendedSystem_box_0 = None

                Coordinates_box_1 = None
                Structure_box_1 = None
                binCoordinates_box_1 = None
                extendedSystem_box_1 = None

        gomc_control.write_gomc_control_file(
            charmm_object_with_files,
            input_name_control_file_name,
            used_ensemble,
            MC_steps,
            temperature,
            ff_psf_pdb_file_directory=None,
            check_input_files_exist=False,
            Parameters="{}.inp".format(ff_filename_str),
            Restart=Restart,
            RestartCheckpoint=True,
            ExpertMode=False,
            Coordinates_box_0=Coordinates_box_0,
            Structure_box_0=Structure_box_0,
            binCoordinates_box_0=binCoordinates_box_0,
            extendedSystem_box_0=extendedSystem_box_0,
            binVelocities_box_0=None,
            Coordinates_box_1=Coordinates_box_1,
            Structure_box_1=Structure_box_1,
            binCoordinates_box_1=binCoordinates_box_1,
            extendedSystem_box_1=extendedSystem_box_1,
            binVelocities_box_1=None,
            input_variables_dict={
                "PRNG": seed_no,
                "Pressure": pressure,
                "Ewald": Ewald,
                "ElectroStatic": ElectroStatic,
                "VDWGeometricSigma": VDWGeometricSigma,
                "Rcut": Rcut,
                "VolFreq": VolFreq[-1],
                "SwapFreq": SwapFreq[-1],
                "DisFreq": DisFreq[-1],
                "RotFreq": RotFreq[-1],
                "RegrowthFreq": RegrowthFreq[-1],
                "OutputName": output_name_control_file_name,
                "EqSteps": EqSteps,
                "PressureCalc": output_true_list_input,
                "RestartFreq": output_true_list_input,
                "CheckpointFreq": output_true_list_input,
                "ConsoleFreq": output_true_list_input,
                "BlockAverageFreq": output_true_list_input,
                "HistogramFreq": output_true_list_input,
                "CoordinatesFreq": output_false_list_input,
                "DCDFreq": output_true_list_input,
                "LRC": LRC,
                "RcutLow": 1,
                "CBMC_First": 12,
                "CBMC_Nth": 10,
                "CBMC_Ang": 50,
                "CBMC_Dih": 50,
            },
        )

        print("#**********************")
        print("Completed: production NPT or GEMC-NVT GOMC control file writing")
        print("#**********************")
        # ******************************************************
        # production NPT or GEMC-NVT - GOMC control file writing  (end)
        # ******************************************************


# ******************************************************
# ******************************************************
# Creating GOMC files (pdb, psf, force field (FF), and gomc control files (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# melt_NVT -starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
@Project.pre(lambda j: j.sp.engine == "gomc")
@Project.pre(part_2a_melt_equilb_NVT_control_file_written)
@Project.post(part_3a_output_melt_equilb_NVT_started)
@Project.post(part_4a_job_melt_equilb_NVT_completed_properly)
@Project.operation.with_directives(
    {
        "np": lambda job: ff_info_dict.get(job.sp.forcefield_name).get(
            "ncpu", 1
        ),
        "ngpu": lambda job: ff_info_dict.get(job.sp.forcefield_name).get(
            "ngpu", 0
        ),
        "memory": memory_needed,
        "walltime": walltime_gomc_hr,
    }
)
@flow.with_job
@flow.cmd
def run_melt_equilb_NVT_gomc_command(job):
    """Run the gomc melt_equilb_NVT simulation."""
    print("#**********************")
    print("# Started the run_melt_NVT_gomc_command.")
    print("#**********************")

    control_file_name_str = melt_equilb_NVT_control_file_name_str

    print(f"Running simulation job id {job}")
    run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
        str(job.doc.gomc_binary_path),
        str(job.doc.melt_NVT_gomc_binary_file),
        str(ff_info_dict.get(job.sp.forcefield_name).get("ncpu")),
        str(control_file_name_str),
        str(control_file_name_str),
    )

    return run_command


# ******************************************************
# ******************************************************
# melt_NVT - including GOMC control file writing and starting the GOMC simulation (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# equilb_NVT - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
@Project.pre(lambda j: j.sp.engine == "gomc")
@Project.pre(part_4a_job_melt_equilb_NVT_completed_properly)
@Project.pre(part_2b_equilb_NVT_control_file_written)
@Project.post(part_3b_output_equilb_NVT_started)
@Project.post(part_4b_job_equilb_NVT_completed_properly)
@Project.operation.with_directives(
    {
        "np": lambda job: ff_info_dict.get(job.sp.forcefield_name).get("ncpu"),
        "ngpu": lambda job: ff_info_dict.get(job.sp.forcefield_name).get(
            "ngpu", 0
        ),
        "memory": memory_needed,
        "walltime": walltime_gomc_hr,
    }
)
@flow.with_job
@flow.cmd
def run_equilb_NVT_gomc_command(job):
    """Run the gomc equilb_NVT simulation."""
    print("#**********************")
    print("# Started the run_NVT_gomc_command.")
    print("#**********************")

    control_file_name_str = equilb_NVT_control_file_name_str

    print(f"Running simulation job id {job}")
    run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
        str(job.doc.gomc_binary_path),
        str(job.doc.equilb_NVT_gomc_binary_file),
        str(ff_info_dict.get(job.sp.forcefield_name).get("ncpu")),
        str(control_file_name_str),
        str(control_file_name_str),
    )

    return run_command


# ******************************************************
# ******************************************************
# equilb_NVT - starting the GOMC simulation (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# Use pymbar to evaluate if the system came to equilibrium based on
# the density (NPT or NVT -> box 0 density, GEMC or GCMC box 1 density) (start)
# ******************************************************
# ******************************************************
def test_pymbar_stabilized_equilb_design_ensemble(job):
    """Test if the simulation has come to equilibrium via pymbar."""
    print("#**********************")
    print("# Started the test_pymbar_stabilized_equilb_design_ensemble")
    print("#**********************")

    fraction_data_required_for_equilbrium = 0.25  # float
    data_points_to_skip_for_equilbrium = 1  # int
    equilb_plot_base_name = "pymbar_equilb_design_ensemble_plot"

    if use_pymbar is True:
        if gomc_sim_completed_properly(
            job,
            job.doc.equilb_design_ensemble_dict[
                str(job.doc.equilb_design_ensemble_number)
            ]["output_name_control_file_name"],
        ):

            box_0_filename = "Blk_{}_BOX_0.dat".format(
                job.doc.equilb_design_ensemble_dict[
                    str(job.doc.equilb_design_ensemble_number)
                ]["output_name_control_file_name"]
            )
            box_0_directory_filename = f"{box_0_filename}"

            read_csv_box_0_data = pd.read_csv(
                box_0_directory_filename,
                sep="\s+",
                header=0,
                na_values="NaN",
                index_col=0,
            )
            read_csv_density_box_0_np_array = np.array(
                read_csv_box_0_data.loc[:, "TOT_DENS"]
            )
            read_csv_total_energy_box_0_np_array = np.array(
                read_csv_box_0_data.loc[:, "TOT_EN"]
            )

            [
                equib_box_0_density_bool,
                t_box_0_density,
                g_box_0_density,
                Neff_box_0_density,
            ] = is_equilibrated(
                read_csv_density_box_0_np_array,
                fraction_data_required_for_equilbrium,
                data_points_to_skip_for_equilbrium,
            )
            [
                equib_box_0_total_energy_bool,
                t_box_0_total_energy,
                g_box_0_total_energy,
                Neff_box_0_total_energy,
            ] = is_equilibrated(
                read_csv_total_energy_box_0_np_array,
                fraction_data_required_for_equilbrium,
                data_points_to_skip_for_equilbrium,
            )

            # needs to be merged before able to use plotting
            plot_data_with_t0_line(
                "{}_box_0_density".format(equilb_plot_base_name),
                read_csv_density_box_0_np_array,
                title="box 0: density (kg/m$^3$) vs. MC cycles * {}"
                "".format(output_data_every_X_MC_cycles),
                threshold=0,
                overwrite=True,
            )
            # needs to be merged before able to use plotting
            plot_data_with_t0_line(
                "{}_box_0_total_energy".format(equilb_plot_base_name),
                read_csv_total_energy_box_0_np_array,
                title="box 0: total potential energy (K) vs. MC cycles * {}"
                "".format(output_data_every_X_MC_cycles),
                threshold=0,
                overwrite=True,
            )

            if job.doc.production_ensemble in ["GEMC_NPT", "GEMC_NVT"]:
                box_1_filename = "Blk_{}_BOX_1.dat".format(
                    job.doc.equilb_design_ensemble_dict[
                        str(job.doc.equilb_design_ensemble_number)
                    ]["output_name_control_file_name"]
                )
                # box_1_directory_filename = f"workspace/{job.id}/{box_1_filename}"
                box_1_directory_filename = f"{box_1_filename}"

                read_csv_box_1_data = pd.read_csv(
                    box_1_directory_filename,
                    sep="\s+",
                    header=0,
                    na_values="NaN",
                    index_col=0,
                )
                read_csv_density_box_1_np_array = np.array(
                    read_csv_box_1_data.loc[:, "TOT_DENS"]
                )
                read_csv_total_energy_box_1_np_array = np.array(
                    read_csv_box_0_data.loc[:, "TOT_EN"]
                )
                read_csv_total_energy_box_0_and_1_np_array = (
                    read_csv_total_energy_box_0_np_array
                    + read_csv_total_energy_box_1_np_array
                )

                [
                    equib_box_1_density_bool,
                    t_box_1_density,
                    g_box_1_density,
                    Neff_box_1_density,
                ] = is_equilibrated(
                    read_csv_density_box_1_np_array,
                    fraction_data_required_for_equilbrium,
                    data_points_to_skip_for_equilbrium,
                )
                [
                    equib_box_1_total_energy_bool,
                    t_box_1_total_energy,
                    g_box_1_total_energy,
                    Neff_box_1_total_energy,
                ] = is_equilibrated(
                    read_csv_total_energy_box_1_np_array,
                    fraction_data_required_for_equilbrium,
                    data_points_to_skip_for_equilbrium,
                )
                [
                    equib_box_0_and_1_total_energy_bool,
                    t_box_0_and_1_total_energy,
                    g_box_0_and_1_total_energy,
                    Neff_box_0_and_1_total_energy,
                ] = is_equilibrated(
                    read_csv_total_energy_box_0_and_1_np_array,
                    fraction_data_required_for_equilbrium,
                    data_points_to_skip_for_equilbrium,
                )
                # needs to be merged before able to use plotting

                plot_data_with_t0_line(
                    "{}_box_1_density".format(equilb_plot_base_name),
                    read_csv_density_box_1_np_array,
                    title="box 1: density (kg/m$^3$) vs. MC cycles * {}"
                    "".format(output_data_every_X_MC_cycles),
                    threshold=0,
                    overwrite=True,
                )
                plot_data_with_t0_line(
                    "{}_box_1_total_energy".format(equilb_plot_base_name),
                    read_csv_total_energy_box_1_np_array,
                    title="box 0: total potential energy (K) vs. MC cycles * {}"
                    "".format(output_data_every_X_MC_cycles),
                    threshold=0,
                    overwrite=True,
                )
                plot_data_with_t0_line(
                    "{}_sum_of_box_0_and_1_total_energy".format(
                        equilb_plot_base_name
                    ),
                    read_csv_total_energy_box_0_and_1_np_array,
                    title="box 0 and 1: summed total potential energy (K) vs. MC cycles * {}"
                    "".format(output_data_every_X_MC_cycles),
                    threshold=0,
                    overwrite=True,
                )

            if job.doc.production_ensemble in ["GCMC", "NVT", "NPT"]:
                if (
                    equib_box_0_density_bool is True
                    and equib_box_0_total_energy_bool is True
                ):
                    job.doc.stable_equilb_design_ensemble = True

            elif job.doc.production_ensemble in ["GEMC_NPT", "GEMC_NVT"]:
                if (
                    equib_box_0_density_bool is True
                    and equib_box_0_total_energy_bool is True
                    and equib_box_1_density_bool is True
                    and equib_box_1_total_energy_bool is True
                    and equib_box_0_and_1_total_energy_bool is True
                ):
                    job.doc.stable_equilb_design_ensemble = True

    elif use_pymbar is False:
        job.doc.stable_equilb_design_ensemble = True

    print("#**********************")
    print("# Completed the test_pymbar_stabilized_equilb_design_ensemble.")
    print("#**********************")


@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@Project.post(part_4c_job_equilb_design_ensemble_completed_properly)
@flow.with_job
def pymbar_stabilized_equilb_design_ensemble(job):
    """Determine if the simulation has come to equilibrium via pymbar and update the doc json file."""
    print("#**********************")
    print(
        "# Running test if the pymbar_stabilized_equilb_design_ensemble is stable."
    )
    print("#**********************")

    return job.doc.stable_equilb_design_ensemble


# ******************************************************
# ******************************************************
# Use pymbar to evaluate if the system came to equilibrium based on
# the density (NPT or NVT -> box 0 density, GEMC or GCMC box 1 density) (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# equilb NPT or GEMC-NVT - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
@Project.pre(lambda j: j.sp.engine == "gomc")
@Project.pre(part_4b_job_equilb_NVT_completed_properly)
@Project.pre(part_2c_equilb_design_ensemble_control_file_written)
@Project.post(part_3c_output_equilb_design_ensemble_started)
@Project.post(part_4c_job_equilb_design_ensemble_completed_properly)
@Project.post(pymbar_stabilized_equilb_design_ensemble)
@Project.operation.with_directives(
    {
        "np": lambda job: ff_info_dict.get(job.sp.forcefield_name).get("ncpu"),
        "ngpu": lambda job: ff_info_dict.get(job.sp.forcefield_name).get(
            "ngpu", 0
        ),
        "memory": memory_needed,
        "walltime": walltime_gomc_hr,
    }
)
@flow.with_job
def run_equilb_ensemble_gomc_command(job):
    """Run the gomc equilb_ensemble simulation."""
    for run_equilb_ensemble_i in range(
        job.doc.equilb_design_ensemble_number, equilb_design_ensemble_max_number
    ):
        print("#**********************")
        print("# Started the run_equilb_ensemble_gomc_command function.")
        print("#**********************")

        if (
            job.doc.equilb_design_ensemble_number
            >= equilb_design_ensemble_max_number
        ):
            job.doc.equilb_design_ensemble_max_number_under_limit = False

        elif (
            job.doc.stable_equilb_design_ensemble is False
            and job.doc.equilb_design_ensemble_max_number_under_limit is True
        ):

            control_file_name_str = job.doc.equilb_design_ensemble_dict[
                str(job.doc.equilb_design_ensemble_number)
            ]["output_name_control_file_name"]

            print(f"Running simulation job id {job}")
            run_command = "{}/{} +p{} {}.conf > out_{}.dat" "".format(
                str(job.doc.gomc_binary_path),
                str(job.doc.equilb_design_ensemble_gomc_binary_file),
                str(ff_info_dict.get(job.sp.forcefield_name).get("ncpu")),
                str(control_file_name_str),
                str(control_file_name_str),
            )

            exec_run_command = subprocess.Popen(
                run_command, shell=True, stderr=subprocess.STDOUT
            )
            os.waitpid(exec_run_command.pid, 0)  # os.WSTOPPED) # 0)

            test_pymbar_stabilized_equilb_design_ensemble(job)

            if job.doc.stable_equilb_design_ensemble is False:
                # need to add equilb_design_ensemble_number by 1 so it is fixed to run the correct job
                # so it is rerun if restarted
                job.doc.equilb_design_ensemble_number += 1


# ******************************************************
# ******************************************************
# equilb NPT or GEMC-NVT - starting the GOMC simulation (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# production run - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
@Project.pre(part_4c_job_equilb_design_ensemble_completed_properly)
@Project.pre(part_2d_production_control_file_written)
@Project.pre(pymbar_stabilized_equilb_design_ensemble)
@Project.post(part_3d_output_production_run_started)
@Project.post(part_4d_job_production_run_completed_properly)
@Project.operation.with_directives(
    {
        "np": lambda job: ff_info_dict.get(job.sp.forcefield_name).get("ncpu"),
        "ngpu": lambda job: ff_info_dict.get(job.sp.forcefield_name).get(
            "ngpu", 0
        ),
        "memory": memory_needed,
        "walltime": walltime_gomc_hr,
    }
)
@flow.with_job
@flow.cmd
def run_production_run_gomc_command(job):
    """Run the gomc production_run simulation."""
    print("#**********************")
    print("# Started the run_production_run_gomc_command function.")
    print("#**********************")

    control_file_name_str = job.doc.production_run_ensemble_dict[
        str(job.doc.equilb_design_ensemble_number)
    ]["input_name_control_file_name"]

    output_file_name_str = job.doc.production_run_ensemble_dict[
        str(job.doc.equilb_design_ensemble_number)
    ]["output_name_control_file_name"]

    print(f"Running simulation job id {job}")
    run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
        str(job.doc.gomc_binary_path),
        str(job.doc.production_ensemble_gomc_binary_file),
        str(ff_info_dict.get(job.sp.forcefield_name).get("ncpu")),
        str(control_file_name_str),
        str(output_file_name_str),
    )

    return run_command


# ******************************************************
# ******************************************************
# production run - starting the GOMC simulation (end)
# ******************************************************
# ******************************************************
