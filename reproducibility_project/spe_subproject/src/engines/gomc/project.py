"""GOMC's single point energy analysis setup for signac, signac-flow, signac-dashboard for this study."""
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
gomc_binary_path = (
    "/Users/brad/Programs/GOMC/GOMC_dev_zero_point_energy_2_28_22/bin"
)

# force field (FF) file for all simulations in that job
# Note: do not add extensions
ff_filename_str = "in_FF"

# initial mosdef structure and coordinates
# Note: do not add extensions
mosdef_structure_box_0_name_str = "mosdef_gomc_zero_point_energy_box_0"
mosdef_structure_box_1_name_str = "mosdef_gomc_zero_point_energy_box_1"

# The production run using the ensemble used for the simulation design, which
# includes the simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
production_control_file_name_str = "production_run"


path_from_job_to_box_inputs = "../../"

walltime_mosdef_hr = 4
walltime_gomc_hr = 4
memory_needed = 16

use_pymbar = True  # True of False

ff_info_dict = {
    "trappe-ua": {
        "ngpu": 0,
        "ncpu": 1,
        "Ewald": False,
        "ElectroStatic": False,
        "VDWGeometricSigma": False,
    },
    "benzene-ua": {
        "ngpu": 0,
        "ncpu": 1,
        "Ewald": False,
        "ElectroStatic": False,
        "VDWGeometricSigma": False,
    },
    "spce": {
        "ngpu": 0,
        "ncpu": 1,
        "Ewald": True,
        "ElectroStatic": True,
        "VDWGeometricSigma": False,
    },
    "oplsaa": {
        "ngpu": 0,
        "ncpu": 1,
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

    job.doc.production_run_ensemble_dict = {
        "input_name_control_file_name": production_control_file_name_str,
        "output_name_control_file_name": production_control_file_name_str,
    }

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

# check if GOMC-MOSDEF wrote the gomc files
# @Project.pre(select_production_ensemble)
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def mosdef_input_written(job):
    """Check that the mosdef files (psf, pdb, and force field (FF) files) are written ."""
    file_written_bool = False

    if job.sp.ensemble in ["NPT", "NVT"]:
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
    elif job.sp.ensemble in ["GCMC", "GEMC_NPT", "GEMC_NPT"]:
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


# checking if the GOMC control file is written for the production run
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_2a_production_control_file_written(job):
    """General check that the production run (run temperature) gomc control file is written."""
    try:
        return gomc_control_file_written(job, production_control_file_name_str)
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


# check if production GOMC run is started by seeing if the GOMC consol file and the merged psf exist
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_3a_output_production_run_started(job):
    """Check to see if the production run (set temperature) gomc simulation is started."""
    try:
        return gomc_simulation_started(job, production_control_file_name_str)

    except:
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


# check if production GOMC run completed by checking the end of the GOMC consol file
@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_4a_job_production_run_completed_properly(job):
    """Check to see if the production run (set temperature) gomc simulation was completed properly."""
    return gomc_sim_completed_properly(job, production_control_file_name_str)


# ******************************************************
# ******************************************************
# check if GOMC simulation are completed properly (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if production energy analysis is completed (start)
# ******************************************************
# ******************************************************


@Project.label
@Project.pre(lambda j: j.sp.engine == "gomc")
@flow.with_job
def part_5a_job_production_energy_analysis_written(job):
    """Check that the production energy analysis is written  ."""
    file_written_bool = False
    # if (job.isfile(f"{path_from_job_to_box_inputs}/log-spe.txt")):
    #    file_written_bool = True

    if job.isfile("log-spe.txt"):
        file_written_bool = True

    return file_written_bool


# ******************************************************
# ******************************************************
# check if production energy analysis is written (end)
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

    if job.sp.molecule in ["waterSPCE", "benzeneUA"]:
        gomc_fix_bonds_angles_list = Molecule_ResName_List
    else:
        gomc_fix_bonds_angles_list = None

    if job.sp.ensemble in ["NVT", "NPT"]:
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

    elif job.sp.ensemble in ["GCMC", "GEMC_NVT", "GEMC_NPT"]:
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
@Project.post(part_2a_production_control_file_written)
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
    # production NPT or GEMC-NVT - GOMC control file writing  (start)
    # ******************************************************

    print("#**********************")
    print("Started: production NPT or GEMC-NVT GOMC control file writing")
    print("#**********************")

    Restart = True

    temperature = (job.sp.temperature * u.K).to_value("K")
    pressure = (job.sp.pressure * u.kPa).to_value("bar")

    input_name_control_file_name = job.doc.production_run_ensemble_dict[
        "input_name_control_file_name"
    ]
    output_name_control_file_name = job.doc.production_run_ensemble_dict[
        "output_name_control_file_name"
    ]

    MC_steps = 2  # set value for paper = 1

    output_data_every_X_MC_steps = 1  # set value for paper = 10
    EqSteps = 1

    seed_no = job.doc.replica_number_int

    production_ensemble = job.sp.ensemble

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
        Potential = "VDW"
        try:
            if job.sp.long_range_correction == "None":
                LRC = False
            elif job.sp.long_range_correction == "energy_pressure":
                LRC = True
            else:
                raise ValueError("ERROR: Not a valid cutoff_style")
        except:
            LRC = False

    elif job.sp.cutoff_style == "shift":
        Potential = "SHIFT"
        LRC = False
    else:
        raise ValueError("ERROR: Not a valid cutoff_style")

    # output all data and calc frequecy
    output_true_list_input = [
        True,
        int(output_data_every_X_MC_steps),
    ]
    output_false_list_input = [
        False,
        int(output_data_every_X_MC_steps),
    ]

    if production_ensemble in ["NVT", "NPT"]:
        used_ensemble = "NPT"

        if job.sp.molecule in ["methaneUA"]:
            VolFreq = (0.01,)
            SwapFreq = (None,)
            DisFreq = (0.99,)
            RotFreq = (None,)
            RegrowthFreq = (None,)

        elif job.sp.molecule in ["pentaneUA", "ethanolAA"]:
            VolFreq = (0.01,)
            SwapFreq = (None,)
            DisFreq = (0.33,)
            RotFreq = (0.33,)
            RegrowthFreq = (0.33,)

        elif job.sp.molecule in ["waterSPCE", "benzeneUA"]:
            VolFreq = (0.01,)
            SwapFreq = (None,)
            DisFreq = (0.49,)
            RotFreq = (0.5,)
            RegrowthFreq = (None,)

        else:
            raise ValueError(
                "Moleules MC move rations not listed in the GOMC control file writer."
            )

        restart_rel_dir = (
            "../../src/system_snapshots/gomc_NPT_percise_coordinates"
        )
        if Restart is True:
            Coordinates_box_0 = "mosdef_gomc_zero_point_energy_box_0.pdb"
            Structure_box_0 = "mosdef_gomc_zero_point_energy_box_0.psf"
            binCoordinates_box_0 = (
                f"{restart_rel_dir}/{job.sp.molecule}"
                f"/mosdef_gomc_zero_point_energy_box_0_restart.coor"
            )
            extendedSystem_box_0 = (
                f"{restart_rel_dir}/{job.sp.molecule}"
                f"/mosdef_gomc_zero_point_energy_box_0_restart.xsc"
            )

        elif Restart is False:
            raise ValueError(
                "ERROR: GOMC zero point energies need to be calculated using a "
                "restart coor and xsc files where were manually created in VMD."
            )

        Coordinates_box_1 = None
        Structure_box_1 = None
        binCoordinates_box_1 = None
        extendedSystem_box_1 = None

    elif job.sp.ensemble in ["GCMC", "GEMC_NVT", "GEMC_NPT"]:
        raise ValueError(
            "ERROR: GOMC zero point energies need to be calculated using and NPT ensemble."
        )

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
            "Potential": Potential,
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
# production run - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
@Project.pre(part_1a_initial_data_input_to_json)
@Project.pre(part_2a_production_control_file_written)
@Project.post(part_3a_output_production_run_started)
@Project.post(part_4a_job_production_run_completed_properly)
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
        "input_name_control_file_name"
    ]

    output_file_name_str = job.doc.production_run_ensemble_dict[
        "output_name_control_file_name"
    ]

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


# ******************************************************
# ******************************************************
# production run - starting the GOMC energy analysis (start)
# ******************************************************
# ******************************************************
@Project.pre(part_4a_job_production_run_completed_properly)
@Project.post(part_5a_job_production_energy_analysis_written)
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
def run_production_energy_analysis(job):
    """Run the production energy analysis function."""
    print("#**********************")
    print("# Started the production energy analysis function.")
    print("#**********************")

    with open("out_production_run.dat", "r") as fp:
        out_gomc = fp.readlines()
        for i, line in enumerate(out_gomc):
            split_line = line.split()
            if len(split_line) == 13:
                if split_line[0] == "ENER_0:" and int(split_line[1]) == 0:
                    # energies in K units
                    total_energy = None
                    potential_energy = float(split_line[2])
                    #  tail_energy = LRC
                    tail_energy = float(split_line[6])
                    total_electrostatics = None
                    # vdw_energy = INTRA(B) + INTRA(NB) + INTER(LJ) +  LRC
                    vdw_energy = (
                        float(split_line[3])
                        + float(split_line[4])
                        + float(split_line[5])
                        + float(split_line[6])
                    )

                    coul_energy = None
                    pair_energy = None
                    bonds_energy = None
                    angles_energy = None
                    dihedrals_energy = None
                    kspace_energy = None

                    intra_bonded = float(split_line[3])
                    intra_non_bonded = float(split_line[4])
                    inter_LJ = float(split_line[5])

    energy_titles = [
        "total_energy",
        "potential_energy",
        "total_electrostatics",
        "vdw_energy",
        "coul_energy",
        "pair_energy",
        "bonds_energy",
        "angles_energy",
        "dihedrals_energy",
        "tail_energy",
        "kspace_energy",
        "intra_bonded",
        "intra_non_bonded",
        "inter_LJ",
    ]

    energy_K = [
        total_energy,
        potential_energy,
        total_electrostatics,
        vdw_energy,
        coul_energy,
        pair_energy,
        bonds_energy,
        angles_energy,
        dihedrals_energy,
        tail_energy,
        kspace_energy,
        intra_bonded,
        intra_non_bonded,
        inter_LJ,
    ]

    energy_kJ_per_mol = []
    for energy_iter in energy_K:
        if isinstance(energy_iter, float):
            energy_value_iter = u.unyt_quantity(float(energy_iter), "K")
            energy_value_iter = energy_value_iter.to_equivalent(
                "kJ/mol", "thermal"
            )
            energy_value_iter = energy_value_iter.to_value("kJ/mol")
            energy_kJ_per_mol.append(energy_value_iter)

        else:
            energy_kJ_per_mol.append(energy_iter)

    energy_kJ_per_mol_df = pd.DataFrame(np.column_stack(energy_kJ_per_mol))
    energy_kJ_per_mol_df.to_csv(
        "log-spe.txt", sep=",", index=False, header=energy_titles
    )


# ******************************************************
# ******************************************************
# production run - starting the GOMC energy analysis (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# signac end code (start)
# ******************************************************
# ******************************************************
if __name__ == "__main__":
    pr = Project()
    pr.main()
# ******************************************************
# ******************************************************
# signac end code (end)
# ******************************************************
# ******************************************************
