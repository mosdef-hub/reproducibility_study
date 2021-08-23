"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib
import sys

import flow
from flow import environments
import numpy as np


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
    return job.isfile("box.lammps")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_copy_files(job):
    return job.isfile("submit.pbs")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_minimized(job):
    return job.isfile("minimized.restart")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_equilibrated_nvt(job):
    return job.isfile("equilibrated_nvt.restart")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_equilibrated_npt(job):
    return job.isfile("equilibrated_npt.restart")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_production(job):
    return job.isfile("production.restart")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_density_data(job):
    return job.isfile("density.dat")


@Project.label
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
def lammps_created_gsd(job):
    return job.isfile("prod.gsd")


# _____________________________________________________________________
"""Setting up workflow operation"""


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.post(lammps_created_box)
@flow.with_job
@flow.cmd
def built_lammps(job):
    # Create a lammps datafile for a specified molecule
    import foyer
    from mbuild.formats.lammpsdata import write_lammpsdata

    sys.path.append(Project().root_directory() + "/..")
    from reproducibility_project.src.molecules.system_builder import (
        construct_system,
    )

    system = construct_system(job.sp)[0]
    parmed_structure = system.to_parmed()
    # Apply forcefield from statepoint
    if job.sp.forcefield_name == "trappe-ua":
        ff = foyer.Forcefield(name="trappe-ua")
    elif job.sp.forcefield_name == "oplsaa":
        ff = foyer.Forcefield(name="oplsaa")
    elif job.sp.forcefield_name == "spce":
        ff = foyer.Forcefield(
            name="spce"
        )  # TODO: Make sure this gets applied correctly
    else:
        raise Exception(
            "No forcefield has been applied to this system {}".format(job.id)
        )
    typed_surface = ff.apply(parmed_structure)
    write_lammpsdata(
        typed_surface,
        "box.lammps",
        atom_style="full",
        unit_style="real",
        mins=[system.get_boundingbox().vectors[0]],
        maxs=[system.get_boundingbox().vectors[1]],
        use_rb_torsions=True,
    )
    return


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(lammps_created_box)
@Project.post(lammps_copy_files)
@flow.with_job
@flow.cmd
def lammps_cp_files(job):
    molecule = job.sp.molecule
    dict_of_lammps_files = {
        "methaneUA": "UAmethane",
        "pentaneUA": "UApentane",
        "benzeneUA": "UAbenzene",
        "waterSPC/E": "SPCEwater",
        "ethanolAA": "AAethanol",
    }

    lmps_submit_path = "../../src/engine_input/lammps/VU_scripts/submit.pbs"
    lmps_run_path = (
        "../../src/engine_input/lammps/input_scripts/in."
        + dict_of_lammps_files[molecule]
    )
    msg = f"cp {lmps_submit_path} {lmps_run_path} ./"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(lammps_copy_files)
@Project.post(lammps_minimized)
@flow.with_job
@flow.cmd
def run_lammps(job):
    dict_of_lammps_files = {
        "methaneUA": "UAmethane",
        "pentaneUA": "UApentane",
        "benzeneUA": "UAbenzene",
        "waterSPC/E": "SPCEwater",
        "ethanolAA": "AAethanol",
    }
    modify_lammps_scripts(
        "in.{}".format(dict_of_lammps_files[job.sp.molecule]), job
    )
    modify_submit_scripts(
        "in.{}".format(dict_of_lammps_files[job.sp.molecule]), job.id[0:3], 8
    )
    msg = f"qsub submit.pbs"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(lammps_production)
@Project.post(lambda j: j.isfile("prod.gsd"))
@flow.with_job
@flow.cmd
def lammps_calc_rdf(job):
    # Create rdf data from the production run
    import mbuild as mb
    import MDAnalysis as mda
    u = mda.Universe("box.lammps", topology_format="DATA")
    u.load_new("prod.xtc")
    u.trajectory[-1]
    parmed_structure = u.atoms.convert_to("PARMED")
    mb.formats.gsdwriter.write_gsd(parmed_structure, "prod.gsd")
    # TODO: Use freud rdf PR to create an RDF from the gsd file
    return

@Project.operation
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(lammps_production)
@Project.post(lambda j: j.isfile("density.csv"))
@flow.with_job
@flow.cmd
def lammps_create_density_csv(job):
    # Create rdf data from the production run
    import lammps_thermo
    thermo = lammps_thermo.load("log.lammps", skip_sections=9)
    density = thermo.prop("Density")
    list_d = np.concatenate(density).ravel() 
    np.savetxt('density.csv', list_d, delimiter=',', fmt='%f')
    return

def modify_submit_scripts(filename, jobid, cores):
    # Modify Submit Scripts
    with open("submit.pbs", "r") as f:
        lines = f.readlines()
        lines[1] = "#PBS -N {}-{}\n".format(filename[3:], jobid)
        lines[11] = "mpirun -np {} lmp < {}\n".format(cores, filename)
    with open("submit.pbs", "w") as f:
        f.writelines(lines)
    return


def modify_lammps_scripts(filename, job):
    with open(filename, "r") as f:
        lines = f.readlines()
        lines[7] = "pair_style     lj/cut/coul/cut {}\n".format(
            job.sp.r_cut * 10
        )  # nm to angstrom
        lines[21] = "variable tsample equal {} #kelvin\n".format(
            job.sp.temperature
        )  # kelvin
        lines[22] = "variable psample equal {} #atm\n".format(
            job.sp.pressure / 101.325
        )  # kPa to atm
        lines[42] = "velocity all create {} {} dist gaussian\n".format(
            job.sp.temperature, job.sp.replica+1
        )
    with open(filename, "w") as f:
        f.writelines(lines)


if __name__ == "__main__":
    pr = Project()
    pr.main()
