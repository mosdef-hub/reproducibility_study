"""Setup for signac, signac-flow, signac-dashboard for running MCCCS-MN simulations for the reproducibility study."""
import fileinput
import math
import os
import pathlib
import shutil
from glob import glob

import flow
from flow import FlowProject, environments
from flow.environment import DefaultSlurmEnvironment


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


class Metropolis(DefaultSlurmEnvironment):  # Grid(StandardEnvironment):
    """Subclass of DefaultSlurmEnvironment for Siepmann Group's Metropolis cluster."""

    hostname_pattern = r".*\.metropolis2\.chem\.umn\.edu"
    template = "metropolis2.sh"

    @classmethod
    def add_args(cls, parser):
        """Add command line arguments to the submit call."""
        parser.add_argument(
            "--walltime",
            type=float,
            default=96,
            help="Walltime for this submission",
        )


ex = Project.make_group(name="ex")


def mc3s_exec():
    """Return the path of MCCCS-MN executable."""
    return "/home/rs/group-code/MCCCS-MN-9-21/exe-ifort-9-21/src/topmon"


def print_running_string(job, step):
    """Print details about the stage that is starting."""
    print(
        "Running {} for {} molecule = {}, ensemble = {}, temperature= {} K, pressure = {} kPa, replica = {}.".format(
            step,
            job,
            job.sp.molecule,
            job.sp.ensemble,
            job.sp.temperature,
            job.sp.pressure,
            job.sp.replica,
        )
    )


def print_completed_string(job, step):
    """Print details about the stage that just completed."""
    print(
        "Completed {} for {} molecule = {}, ensemble = {}, temperature= {} K, pressure = {} kPa, replica = {}.".format(
            step,
            job,
            job.sp.molecule,
            job.sp.ensemble,
            job.sp.temperature,
            job.sp.pressure,
            job.sp.replica,
        )
    )


"""Setting progress label"""


@Project.label
def has_fort_files(job):
    """Check if the job has all four equired fort.4 files."""
    return (
        job.isfile("fort.4.melt")
        and job.isfile("fort.4.cool")
        and job.isfile("fort.4.equil")
        and job.isfile("fort.4.prod")
    )


@Project.label
@Project.pre(lambda j: j.sp.engine == "mcccs")
def files_ready(job):
    """Check if the keywords in the fort.4 files have been replaced."""
    # Link: https://stackoverflow.com/questions/32749350/check-if-a-string-is-in-a-file-with-python
    if job.sp.ensemble == "GEMC-NVT":
        job.doc.files_ready = False
        file_names = ["melt", "cool", "equil", "prod"]
        keywords = [
            "NCHAIN1",
            "NCHAIN2",
            "LENGTH1",
            "LENGTH2",
            "TEMPERATURE",
            "RCUT",
            "NCHAINTOT",
            "INIX1",
            "INIY1",
            "INIZ1",
            "INIX2",
            "INIY2",
            "INIZ2",
            "VARIABLES",
        ]
        c = 0
        for name in file_names:
            file_name = job.ws + "/fort.4." + name
            i = 0
            if not job.isfile("fort.4." + name):
                continue
            for i in range(len(keywords)):
                with open(file_name) as myfile:
                    if keywords[i] in myfile.read():
                        c += 1
        if c == 0:
            job.doc.files_ready = True

        return job.doc.files_ready

    if job.sp.ensemble == "NPT":
        job.doc.files_ready = False
        file_names = ["melt", "cool", "equil", "prod"]
        keywords = [
            "NCHAIN",
            "LENGTH",
            "TEMPERATURE",
            "PRESSURE",
            "RCUT",
            "VARIABLES",
        ]
        c = 0
        for name in file_names:
            file_name = job.ws + "/fort.4." + name
            i = 0
            if not job.isfile("fort.4." + name):
                continue
            for i in range(len(keywords)):
                with open(file_name) as myfile:
                    if keywords[i] in myfile.read():
                        c += 1
        if c == 0:
            job.doc.files_ready = True

        return job.doc.files_ready


@Project.label
def has_restart_file(job):
    """Check if the job has a restart file."""
    return job.isfile("fort.77")


@Project.label
def has_topmon(job):
    """Check if the job has a topmon (FF) file."""
    return job.isfile("topmon.inp")


@Project.label
def equil_replicate_set(job):
    """Check if number of equil replicates done has been set."""
    return isinstance(job.doc.get("equil_replicates_done"), int)


@Project.label
def equil_replicate_set(job):
    """Check if number of equil replicates done has been set."""
    try:
        return isinstance(job.doc.equil_replicates_done, int)
    except AttributeError:
        return False


@Project.label
def replicate_set(job):
    """Check if number of replicates for prod has been set."""
    return (job.doc.get("num_prod_replicates") == 4) and isinstance(
        job.doc.get("prod_replicates_done"), int
    )


@Project.label
def all_prod_replicates_done(job):
    """Check if all prod replicate simulations completed."""
    try:
        a = job.doc.get("prod_replicates_done")
        b = job.doc.get("num_prod_replicates")
        if a >= b:
            print("simulation complete for {} job".format(job))
        return a >= b
    except (AttributeError, KeyError) as e:
        return False


@Project.label
def melt_finished(job):
    """Check if melt stage is finished."""
    step = "melt"
    run_file = job.ws + "/run.{}".format(step)
    if job.isfile("run.{}".format(step)):
        with open(run_file) as myfile:
            return "Program ended" in myfile.read()
    else:
        return False


@Project.label
def cool_finished(job):
    """Check if cool stage is finished."""
    step = "cool"
    run_file = job.ws + "/run.{}".format(step)
    if job.isfile("run.{}".format(step)):
        with open(run_file) as myfile:
            return "Program ended" in myfile.read()
    else:
        return False


@Project.label
def equil_finished(job):
    """Check if equil stage is finished."""
    try:

        step = "equil" + str(job.doc.equil_replicates_done - 1)
    except (KeyError, AttributeError):
        step = "equil" + "0"

    run_file = job.ws + "/run.{}".format(step)
    if job.isfile("run.{}".format(step)):
        with open(run_file) as myfile:
            return "Program ended" in myfile.read()
    else:
        return False


def sanitize_npt_log(step):
    """Sanitize the output logs for NPT simulations."""
    import numpy as np

    files = sorted(glob("fort*12*{}*".format(step)))
    arrays = []
    for filecurrent in files:
        array = np.genfromtxt(filecurrent, skip_header=1)
        arrays.append(array)
    arrays = np.vstack(arrays)
    arrays[:, 0] = arrays[:, 0] / 10  # Ang to nm
    arrays[:, 1] = arrays[:, 1] / 10  # Ang to nm
    arrays[:, 2] = arrays[:, 2] / 10  # Ang to nm
    arrays[:, 3] = arrays[:, 3] * 0.008314410016255453  # K to kJ/mol
    arrays[:, 4] = arrays[:, 4]  # Pressure kPa to kPa
    density_array = arrays[:, 5] / (arrays[:, 0]) ** 3  # density(molecules/nm3)
    density_array = density_array.reshape(density_array.shape[0], 1)
    arrays = np.append(arrays, density_array, axis=1)  # density molcules/nm3

    np.savetxt(
        "{}_log.txt".format(step),
        arrays,
        header="a (nm) \t b (nm) \t c (nm) \t Energy (kJ/mol) \t Pressure (kPa) \t #molecules \t density (molecules/nm3)",
    )
    return arrays


def sanitize_gemc_log(step):
    """Sanitize the output logs for gemc simulations."""
    import numpy as np

    files = sorted(glob("fort*12*{}*".format(step)))
    arrays_box1 = []
    arrays_box2 = []

    for filecurrent in files:
        array = np.genfromtxt(filecurrent, skip_header=1)
        array1 = array[::2, :]
        array2 = array[1::2, :]
        arrays_box1.append(array1)
        arrays_box2.append(array2)
    arrays_box1 = np.vstack(arrays_box1)
    arrays_box2 = np.vstack(arrays_box2)
    arrays_box1[:, 0] = arrays_box1[:, 0] / 10  # Ang to nm
    arrays_box1[:, 1] = arrays_box1[:, 1] / 10  # Ang to nm
    arrays_box1[:, 2] = arrays_box1[:, 2] / 10  # Ang to nm
    arrays_box1[:, 3] = arrays_box1[:, 3] * 0.008314410016255453  # K to kJ/mol
    arrays_box1[:, 4] = arrays_box1[:, 4]  # Pressure kPa to kPa
    density_array = arrays_box1[:, 5] / (arrays_box1[:, 0]) ** 3
    density_array = density_array.reshape(density_array.shape[0], 1)
    arrays_box1 = np.append(
        arrays_box1, density_array, axis=1
    )  # density molcules/nm3
    # arrays_box1[:, 6] = arrays_box1[:, 5]/ (arrays_box1[:, 0]) **3 # density molcules/nm3
    arrays_box2[:, 0] = arrays_box2[:, 0] / 10  # Ang to nm
    arrays_box2[:, 1] = arrays_box2[:, 1] / 10  # Ang to nm
    arrays_box2[:, 2] = arrays_box2[:, 2] / 10  # Ang to nm
    arrays_box2[:, 3] = arrays_box2[:, 3] * 0.008314410016255453  # K to kJ/mol
    arrays_box2[:, 4] = arrays_box2[:, 4]  # Pressure kPa to kPa
    # arrays_box2[:, 6] = arrays_box2[:, 5]/ (arrays_box2[:, 0]) **3 # density molcules/nm3
    density_array = arrays_box2[:, 5] / (arrays_box2[:, 0]) ** 3
    density_array = density_array.reshape(density_array.shape[0], 1)
    arrays_box2 = np.append(
        arrays_box2, density_array, axis=1
    )  # density molcules/nm3

    np.savetxt(
        "{}_log_box1.txt".format(step),
        arrays_box1,
        header="a (nm) \t b (nm) \t c (nm) \t Energy (kJ/mol) \t Pressure (kPa) \t #molecules \t density (molecules/nm^3)",
    )
    np.savetxt(
        "{}_log_box2.txt".format(step),
        arrays_box2,
        header="a (nm) \t b (nm) \t c (nm) \t Energy (kJ/mol) \t Pressure (kPa) \t #molecules \t density (molecules/nm^3)",
    )
    return arrays_box1, arrays_box2


def system_equilibrated(job):
    """Check if the system is equilibrated."""
    from reproducibility_project.src.analysis.equilibration import (
        is_equilibrated,
    )

    # If a system is already equilibrated, we don't want to check equilibration again, so read information from the file and return True if equilibrated.
    if job.doc.get("is_equilibrated") == True:
        with open(job.ws + "/equil_information.txt", "r") as f:
            print(f.read())
        return True

    with job:
        files = glob("fort*12*{}*".format("equil"))

        if len(files) < 2:  # at least do two loops of equilibration

            print(
                "equils done is less than 2 for {} molecule = {}, ensemble = {}, temperature= {} K, pressure = {} kPa.".format(
                    job,
                    job.sp.molecule,
                    job.sp.ensemble,
                    job.sp.temperature,
                    job.sp.pressure,
                )
            )
            return False

        if job.sp.ensemble == "NPT":
            equil_log = sanitize_npt_log("equil")
            # Now run pymbar on box length and box energy
            equil_status_density = is_equilibrated(
                equil_log[:, 6], threshold_fraction=0.2, nskip=100
            )
            equil_status_energy = is_equilibrated(
                equil_log[:, 3], threshold_fraction=0.2, nskip=100
            )
            if (equil_status_density[0] and equil_status_energy[0]) == False:
                print(
                    "System {} is not equilibrated. Completed {} equil loops".format(
                        job, job.doc.get("equil_replicates_done")
                    )
                )
                print("Equil status density1={}".format(equil_status_density))
                print("Equil status energy1={}".format(equil_status_energy))
                if len(files) >= 3:
                    print(
                        "Even though the system is not equilibrated according to pymbar, we are considering this system equilibrated as 3 equil loops are completed"
                    )
                    text_file = open("equil_information.txt", "w")
                    n = text_file.write(
                        "Even though the system is not equilibrated according to pymbar, we are considering this system equilibrated as 3 equil loops are completed"
                    )
                    text_file.close()

                    job.doc.is_equilibrated = True
                    return True
                return False
            if (equil_status_density[0] and equil_status_energy[0]) == True:
                print(
                    "System {} is equilibrated at cycle {}. Completed {} equil loops".format(
                        job,
                        max(equil_status_density[1], equil_status_energy[1]),
                        job.doc.get("equil_replicates_done"),
                    )
                )
                text_file = open("equil_information.txt", "w")
                n = text_file.write(
                    "System {} is equilibrated at cycle {}. Completed {} equil loops".format(
                        job,
                        max(equil_status_density[1], equil_status_energy[1]),
                        job.doc.get("equil_replicates_done"),
                    )
                )
                text_file.close()
                job.doc.is_equilibrated = True
                return True

        if job.sp.ensemble == "GEMC-NVT":
            print("Checking eqlb for GEMC-NVT")
            equil_log_box1 = sanitize_gemc_log("equil")[0]
            equil_log_box2 = sanitize_gemc_log("equil")[1]
            equil_status_density1 = is_equilibrated(
                equil_log_box1[:, 6], threshold_fraction=0.2, nskip=100
            )
            equil_status_energy1 = is_equilibrated(
                equil_log_box1[:, 3], threshold_fraction=0.2, nskip=100
            )
            equil_status_density2 = is_equilibrated(
                equil_log_box2[:, 6], threshold_fraction=0.2, nskip=100
            )
            equil_status_energy2 = is_equilibrated(
                equil_log_box2[:, 3], threshold_fraction=0.2, nskip=100
            )
            equil_status_total_energy = is_equilibrated(
                equil_log_box1[:, 3] + equil_log_box2[:, 3],
                threshold_fraction=0.2,
                nskip=100,
            )
            if (
                equil_status_density1[0]
                and equil_status_energy1[0]
                and equil_status_density2[0]
                and equil_status_energy2[0]
                and equil_status_total_energy[0]
            ) == False:
                print(
                    "System {} is not equilibrated. Completed {} equil loops".format(
                        job, job.doc.get("equil_replicates_done")
                    )
                )
                print("Equil status density1={}".format(equil_status_density1))
                print("Equil status energy1={}".format(equil_status_energy1))
                print("Equil status density2={}".format(equil_status_density2))
                print("Equil status energy2={}".format(equil_status_energy2))
                print(
                    "Equil status total energy={}".format(
                        equil_status_total_energy
                    )
                )
                if len(files) >= 3:
                    print(
                        "Even though the system is not equilibrated according to pymbar, we are considering this system equilibrated as 3 equil loops are completed"
                    )
                    text_file = open("equil_information.txt", "w")
                    n = text_file.write(
                        "Even though the system is not equilibrated according to pymbar, we are considering this system equilibrated as 3 equil loops are completed"
                    )
                    text_file.close()
                    job.doc.is_equilibrated = True
                    return True

                return False
            if (
                equil_status_density1[0]
                and equil_status_energy1[0]
                and equil_status_density2[0]
                and equil_status_energy2[0]
                and equil_status_total_energy[0]
            ) == True:
                print(
                    "System {} is equilibrated at cycle {}. Completed {} equil loops".format(
                        job,
                        max(
                            equil_status_density1[1],
                            equil_status_energy1[1],
                            equil_status_density2[1],
                            equil_status_energy2[1],
                            equil_status_total_energy[1],
                        ),
                        job.doc.get("equil_replicates_done"),
                    )
                )
                text_file = open("equil_information.txt", "w")
                n = text_file.write(
                    "System {} is equilibrated at cycle {}. Completed {} equil loops".format(
                        job,
                        max(
                            equil_status_density1[1],
                            equil_status_energy1[1],
                            equil_status_density2[1],
                            equil_status_energy2[1],
                            equil_status_total_energy[1],
                        ),
                        job.doc.get("equil_replicates_done"),
                    )
                )
                text_file.close()
                return True

        if job.sp.ensemble == "GEMC-NVT":
            equil_log_box1 = sanitize_gemc_log("equil")[0]
            equil_log_box2 = sanitize_gemc_log("equil")[1]
            equil_status_length1 = is_equilibrated(
                equil_log_box1[:, 0], threshold=0.2, nskip=100
            )
            equil_status_energy1 = is_equilibrated(
                equil_log_box1[:, 3], threshold=0.2, nskip=100
            )
            equil_status_length2 = is_equilibrated(
                equil_log_box2[:, 0], threshold=0.2, nskip=100
            )
            equil_status_energy2 = is_equilibrated(
                equil_log_box2[:, 3], threshold=0.2, nskip=100
            )
            equil_status_total_energy = is_equilibrated(
                equil_log_box1[:, 3] + equil_log_box2[:, 3],
                threshold=0.2,
                nskip=100,
            )
            if (
                equil_status_length1[0]
                and equil_status_energy1[0]
                and equil_status_length2[0]
                and equil_status_energy2[0]
                and equil_status_total_energy
            ) == False:
                print(
                    "System {} is not equilibrated. Completed {} equil loops".format(
                        job, job.doc.get("equil_replicates_done")
                    )
                )
                return False
            if (
                equil_status_length1[0]
                and equil_status_energy1[0]
                and equil_status_length2[0]
                and equil_status_energy2[0]
                and equil_status_total_energy
            ) == True:
                print(
                    "System {} is equilibrated at cycle {}. Completed {} equil loops".format(
                        job,
                        max(
                            equil_status_length1[1],
                            equil_status_energy1[1],
                            equil_status_length2[1],
                            equil_status_energy2[1],
                            equil_status_total_energy[1],
                        ),
                        job.doc.get("equil_replicates_done"),
                    )
                )
                text_file.close()
                job.doc.is_equilibrated = True
                return True


@Project.label
def prod_finished(job):
    """Check if prod stage is finished."""
    try:

        step = "prod" + str(job.doc.prod_replicates_done - 1)
    except (KeyError, AttributeError):
        step = "prod" + "0"
    run_file = job.ws + "/run.{}".format(step)
    if job.isfile("run.{}".format(step)):
        with open(run_file) as myfile:
            return "Program ended" in myfile.read()
    else:
        return False


"""Setting up workflow operation"""


@Project.operation
@Project.pre(lambda j: j.sp.engine == "mcccs")
@Project.post(
    lambda j: (
        j.isfile("init1.pdb")
        and j.isfile("init1.mol2")
        and j.sp.ensemble == "NPT"
    )
    or (
        (
            j.isfile("init1.pdb")
            and j.isfile("init2.pdb")
            and j.isfile("init1.mol2")
            and j.isfile("init2.mol2")
            and j.sp.ensemble == "GEMC-NVT"
        )
    )
)
@flow.with_job
def save_top(job):
    """Save topology files for the two boxes."""
    from reproducibility_project.src.molecules.system_builder import (
        construct_system,
    )
    from reproducibility_project.src.utils.forcefields import load_ff

    # Create a Compound and save to pdb
    system = construct_system(job.sp)
    ff = load_ff(job.sp.forcefield_name)
    param_system = ff.apply(system[0])
    param_system.save(
        "init1.pdb",
        overwrite=True,
    )
    param_system.save(
        "init1.mol2",
        overwrite=True,
    )

    if system[1] is None:
        return
    param_system = ff.apply(system[1])
    param_system.save(
        "init2.pdb",
        overwrite=True,
    )
    param_system.save(
        "init2.mol2",
        overwrite=True,
    )


@Project.operation
@Project.pre(lambda j: j.sp.engine == "mcccs")
@Project.post(equil_replicate_set)
def set_equil_replicates(job):
    """Copy the files for simulation from engine_input folder."""
    print("equil replicates set for job {}".format(job))
    job.doc.equil_replicates_done = 0


@Project.operation
@Project.pre(lambda j: j.sp.engine == "mcccs")
@Project.post(replicate_set)
def set_prod_replicates(job):
    """Copy the files for simulation from engine_input folder."""
    print("prod replicates set for job {}".format(job))
    job.doc.num_prod_replicates = 4
    job.doc.prod_replicates_done = 0


@ex
@Project.operation
@Project.pre(lambda j: j.sp.engine == "mcccs")
@Project.post(has_fort_files)
def copy_files(job):
    """Copy the files for simulation from engine_input folder."""
    for file in glob(
        Project().root_directory()
        + "/src/engine_input/mcccs/{}/{}/fort.4.*".format(
            job.sp.molecule, job.sp.ensemble
        )
    ):
        shutil.copy(file, job.workspace() + "/")


@ex
@Project.operation
@Project.pre(lambda j: j.sp.engine == "mcccs")
@Project.post(has_topmon)
def copy_topmon(job):
    """Copy topmon.inp from root directory to mcccs directory."""
    shutil.copy(
        Project().root_directory()
        + "/src/engine_input/mcccs/{}/{}/topmon.inp".format(
            job.sp.molecule, job.sp.ensemble
        ),
        job.workspace() + "/",
    )


@ex
@Project.operation
@Project.pre(lambda j: j.sp.engine == "mcccs" and j.sp.ensemble == "NPT")
@Project.pre(has_fort_files)
@Project.post(files_ready)
def replace_keyword_fort_files_npt(job):
    """Replace keywords with the values of the variables defined in signac statepoint."""
    file_names = ["melt", "cool", "equil", "prod"]
    seed = job.sp.replica
    nchain = job.sp.N_liquid
    length = job.sp.box_L_liq * 10  # nm to A
    temperature = job.sp.temperature
    pressure = job.sp.pressure / 1000  # kPa to MPa
    rcut = job.sp.r_cut * 10
    variables = [nchain, length, temperature, pressure, seed, rcut]
    keywords = ["NCHAIN", "LENGTH", "TEMPERATURE", "PRESSURE", "SEED", "RCUT"]
    for name in file_names:
        file_name = job.ws + "/fort.4." + name
        i = 0
        for i in range(len(variables)):
            with fileinput.FileInput(file_name, inplace=True) as file:
                for line in file:
                    print(line.replace(keywords[i], str(variables[i])), end="")


@ex
@Project.operation
@Project.pre(lambda j: j.sp.engine == "mcccs" and j.sp.ensemble == "GEMC-NVT")
@Project.pre(has_fort_files)
@Project.post(files_ready)
def replace_keyword_fort_files_gemc(job):
    """Replace keywords with the values of the variables defined in signac statepoint."""
    file_names = ["melt", "cool", "equil", "prod"]
    seed = job.sp.replica
    nchain1 = job.sp.N_liquid
    length1 = job.sp.box_L_liq * 10  # nm to A
    nchain2 = job.sp.N_vap
    length2 = job.sp.box_L_vap * 10  # nm to A
    temperature = job.sp.temperature
    pressure = job.sp.pressure / 1000  # kPa to MPa
    rcut = job.sp.r_cut * 10
    nchaintot = nchain1 + nchain2
    inix1 = 1 + math.ceil(nchain1 ** 0.33)
    iniy1 = 1 + math.ceil(nchain1 ** 0.33)
    iniz1 = 1 + math.ceil(nchain1 ** 0.33)
    inix2 = 1 + math.ceil(nchain2 ** 0.33)
    iniy2 = 1 + math.ceil(nchain2 ** 0.33)
    iniz2 = 1 + math.ceil(nchain2 ** 0.33)
    variables = [
        nchain1,
        length1,
        nchain2,
        length2,
        temperature,
        pressure,
        seed,
        rcut,
        nchaintot,
        inix1,
        iniy1,
        iniz1,
        inix2,
        iniy2,
        iniz2,
    ]
    keywords = [
        "NCHAIN1",
        "LENGTH1",
        "NCHAIN2",
        "LENGTH2",
        "TEMPERATURE",
        "PRESSURE",
        "SEED",
        "RCUT",
        "NCHAINTOT",
        "INIX1",
        "INIY1",
        "INIZ1",
        "INIX2",
        "INIY2",
        "INIZ2",
    ]
    for name in file_names:
        file_name = job.ws + "/fort.4." + name
        i = 0
        for i in range(len(variables)):
            with fileinput.FileInput(file_name, inplace=True) as file:
                for line in file:
                    print(line.replace(keywords[i], str(variables[i])), end="")


@ex
@Project.operation
@Project.pre(lambda j: j.sp.engine == "mcccs")
@Project.post(has_restart_file)
def make_restart_file(job):
    """Make a restart file for the job using fort77maker."""
    from reproducibility_project.src.molecules.system_builder import (
        construct_system,
        get_molecule,
    )

    molecules = [get_molecule(job.sp)]
    system = construct_system(job.sp, constrain=True)
    """Make a fort77 file for the job."""
    if job.sp.ensemble == "GEMC-NVT":
        #        from fort77maker_twobox import fort77writer
        from reproducibility_project.src.engine_input.mcccs.fort77maker_twobox import (
            fort77writer,
        )

        filled_boxes = system
        fort77writer(
            molecules,
            filled_boxes,
            output_file=job.ws + "/fort.77",
            xyz_file=[
                job.ws + "/initial_structure_box1.xyz",
                job.ws + "/initial_structure_box2.xyz",
            ],
        )

        return
    elif job.sp.ensemble == "NPT":
        # from fort77maker_onebox import fort77writer
        from reproducibility_project.src.engine_input.mcccs.fort77maker_onebox import (
            fort77writer,
        )

        filled_box = system[0]
        print("The filled box in make_restart file is {}".format(filled_box))
        fort77writer(
            molecules,
            filled_box,
            output_file=job.ws + "/fort.77",
            xyz_file=job.ws + "/initial_structure.xyz",
        )
        return


@ex
@Project.operation
@Project.pre(lambda j: j.sp.engine == "mcccs")
@Project.pre(has_restart_file)
@Project.pre(has_fort_files)
@Project.pre(has_topmon)
@Project.pre(files_ready)
@Project.post(melt_finished)
def run_melt(job):
    """Run melting stage."""
    from subprocess import PIPE, Popen

    step = "melt"
    """Run the melting stage."""
    print_running_string(job, step)
    execommand = mc3s_exec()
    with job:
        shutil.copyfile("fort.4.{}".format(step), "fort.4")
        process = Popen(
            execommand,
            shell=True,
            universal_newlines=True,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE,
        )
        output, error = process.communicate()
        print(output)
        shutil.move("fort.12", "fort.12.{}".format(step))
        shutil.move("box1config1a.xyz", "box1config1a.xyz.{}".format(step))
        shutil.move("run1a.dat", "run.{}".format(step))
        shutil.copy("config1a.dat", "fort.77")
        shutil.move("config1a.dat", "config1a.dat.{}".format(step))
        shutil.move("box1movie1a.pdb", "box1movie1a.{}.pdb".format(step))
        shutil.move("box1movie1a.xyz", "box1movie1a.{}.xyz".format(step))
        print_completed_string(job, step)


@Project.operation
@Project.pre(lambda j: j.sp.engine == "mcccs")
@Project.pre(has_restart_file)
@Project.pre(melt_finished)
@Project.post(cool_finished)
def run_cool(job):
    """Run cool stage."""
    from subprocess import PIPE, Popen

    step = "cool"
    """Run the melting stage."""
    print_running_string(job, step)
    execommand = mc3s_exec()
    with job:

        shutil.copyfile("fort.4.{}".format(step), "fort.4")
        process = Popen(
            execommand,
            shell=True,
            universal_newlines=True,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE,
        )
        output, error = process.communicate()
        print(output)
        shutil.move("fort.12", "fort.12.{}".format(step))
        shutil.move("box1config1a.xyz", "box1config1a.xyz.{}".format(step))
        shutil.move("run1a.dat", "run.{}".format(step))
        shutil.copy("config1a.dat", "fort.77")
        shutil.move("config1a.dat", "config1a.dat.{}".format(step))
        shutil.move("box1movie1a.pdb", "box1movie1a.{}.pdb".format(step))
        shutil.move("box1movie1a.xyz", "box1movie1a.{}.xyz".format(step))
    print_completed_string(job, step)


@Project.operation.with_directives({"walltime": 200})
@Project.pre(lambda j: j.sp.engine == "mcccs")
@Project.pre(has_restart_file)
@Project.pre(cool_finished)
@Project.post(equil_finished)
@Project.post(system_equilibrated)
def run_equil(job):
    """Run equilibration."""
    from subprocess import PIPE, Popen

    step = "equil" + str(job.doc.equil_replicates_done)
    """Run the  equil  stage."""
    print_running_string(job, step)
    execommand = mc3s_exec()
    with job:
        shutil.copyfile("fort.4.equil", "fort.4")
        process = Popen(
            execommand,
            shell=True,
            universal_newlines=True,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE,
        )
        output, error = process.communicate()
        print(output)
        shutil.move("fort.12", "fort.12.{}".format(step))
        shutil.move("box1config1a.xyz", "box1config1a.xyz.{}".format(step))
        shutil.move("run1a.dat", "run.{}".format(step))
        shutil.copy("config1a.dat", "fort.77")
        shutil.move("config1a.dat", "config1a.dat.{}".format(step))
        shutil.move("box1movie1a.pdb", "box1movie1a.{}.pdb".format(step))
        shutil.move("box1movie1a.xyz", "box1movie1a.{}.xyz".format(step))
        if job.sp.ensemble == "GEMC-NVT":
            shutil.move("box2movie1a.pdb", "box1movie1a.{}.pdb".format(step))
            shutil.move("box2movie1a.xyz", "box1movie1a.{}.xyz".format(step))
        job.doc.equil_replicates_done += 1
    print_completed_string(job, step)


@Project.operation.with_directives({"walltime": 200})
@Project.pre(lambda j: j.sp.engine == "mcccs")
@Project.pre(has_restart_file)
@Project.pre(system_equilibrated)
@Project.post(prod_finished)
@Project.post(all_prod_replicates_done)
def run_prod(job):
    """Run production."""
    from subprocess import PIPE, Popen

    replicate = job.doc.prod_replicates_done
    step = "prod" + str(replicate)
    """Run the prod stage."""
    print_running_string(job, step)
    execommand = mc3s_exec()
    with job:
        shutil.copyfile("fort.4.prod", "fort.4")
        process = Popen(
            execommand,
            shell=True,
            universal_newlines=True,
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE,
        )
        output, error = process.communicate()
        print(output)
        shutil.move("fort.12", "fort.12.{}".format(step))
        shutil.move("box1config1a.xyz", "box1config1a.xyz.{}".format(step))
        shutil.move("run1a.dat", "run.{}".format(step))
        shutil.copy("config1a.dat", "fort.77")
        shutil.move("config1a.dat", "config1a.dat.{}".format(step))
        shutil.move("box1movie1a.pdb", "box1movie1a.{}.pdb".format(step))
        shutil.move("box1movie1a.xyz", "box1movie1a.{}.xyz".format(step))
        if job.sp.ensemble == "GEMC-NVT":
            shutil.move("box2movie1a.pdb", "box1movie1a.{}.pdb".format(step))
            shutil.move("box2movie1a.xyz", "box1movie1a.{}.xyz".format(step))

        job.doc.prod_replicates_done += 1
        print_completed_string(job, step)
        if all_prod_replicates_done(job):
            print(
                "All prod replicates done. Simulation finished for {} molecule = {}, ensemble = {}, temperature= {} K, pressure = {} kPa.".format(
                    job,
                    job.sp.molecule,
                    job.sp.ensemble,
                    job.sp.temperature,
                    job.sp.pressure,
                )
            )
            text_file = open("production_information.txt", "w")
            prod_file = text_file.write(
                "All prod replicates done. Simulation finished for {} molecule = {}, ensemble = {}, temperature= {} K, pressure = {} kPa.".format(
                    job,
                    job.sp.molecule,
                    job.sp.ensemble,
                    job.sp.temperature,
                    job.sp.pressure,
                )
            )
            text_file.close()


if __name__ == "__main__":
    pr = Project()
    pr.main()

