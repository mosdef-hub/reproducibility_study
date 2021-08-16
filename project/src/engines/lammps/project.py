"""Setup for signac, signac-flow, signac-dashboard for this study."""
# import foyer
import os
import pathlib

import flow
from flow import environments


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"

#____________________________________________________________________________
"""Setting progress label"""
@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
def lammps_created_box(job):
    return job.isfile("box.lammps")

@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
def lammps_copy_files(job):
    return job.isfile("submit.pbs")

@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
def lammps_minimized(job):
    return job.isfile("minimized.restart")

@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
def lammps_equilibrated_nvt(job):
    return job.isfile("equilibrated_nvt.restart")

@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
def lammps_equilibrated_npt(job):
    return job.isfile("equilibrated_npt.restart")

@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
def lammps_production(job):
    return job.isfile("production.restart")

@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
def lammps_density_data(job):
    return job.isfile("density.dat")

@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
def lammps_created_gsd(job):
    return job.isfile("prod.gsd")

#_____________________________________________________________________
"""Setting up workflow operation"""

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
@Project.post(lammps_created_box)
@flow.cmd
@flow.with_job
def built_lammps(job):
    # Create a lammps datafile for a specified molecule
    from mbuild.formats.lammpsdata import write_lammpsdata

    from project.src.molecules.system_builder import SystemBuilder
    system = SystemBuilder(job)
    parmed_structure = system.to_parmed()
    trappe = foyer.Forcefield(name="trappe-ua")
    typed_surface = trappe.apply(parmed_structure)
    write_lammpsdata(system, "box.lammps", atom_style="full", unit_style="lj", detect_forcefield_style=True)
    return

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
@Project.pre(lammps_created_box)
@Project.post(lammps_copy_files)
@flow.cmd
@flow.with_job
def lammps_cp_files(job):
    molecule = job.sp.molecule
    dict_of_lammps_files = {
                            "methaneUA": "UAmethane"
                            "pentaneUA": "UApentane"
                            "benzeneUA": "UAbenzene"
                            "waterSPC/E": "SPCEwater"
                            "ethanolAA": "AAethanol"
                           }

    lmps_submit_path = "../../engine_input/lammps/submission_scripts/submit.pbs"
    lmps_run_path = "../../engine_input/lammps/submission_scripts/in." + dict_of_lammps_files[molecule]
    msg = f"cp {lmps_inpt_path} {lmps_run_path} ./"
    return msg

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
@Project.pre(lammps_copy_files)
@Project.post(lammps_minimized)
@flow.cmd
@flow.with_job
def lammps_em(job):
    modify_submit_scripts("in.em", str(job.sp.molecule), 8)
    msg = f"qsub submit.pbs"
    return msg

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
@Project.pre(lammps_minimized)
@Project.post(lammps_equilibrated_nvt)
@flow.cmd
@flow.with_job
def lammps_nvt(job):
    modify_submit_scripts("in.nvt", str(job.sp.molecule), 8)
    msg = f"qsub submit.pbs"
    return msg

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
@Project.pre(lammps_equilibrated_nvt)
@Project.post(lammps_equilibrated_npt)
@flow.cmd
@flow.with_job
def lammps_npt(job):
    modify_submit_scripts("in.npt", str(job.sp.molecule), 8)
    msg = f"qsub submit.pbs"
    return msg

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
@Project.pre(lammps_equilibrated_npt)
@Project.post(lammps_production)
@flow.cmd
@flow.with_job
def lammps_prod(job):
    modify_submit_scripts("in.prod", str(job.sp.molecule), 8)
    msg = f"qsub submit.pbs"
    return msg

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
@Project.pre(lammps_equilibrated_production)
@flow.cmd
@flow.with_job
def lammps_calc_density(job):
    # Create a density datafile from the production run
    return

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps")
@Project.pre(lammps_equilibrated_production)
@flow.cmd
@flow.with_job
def lammps_calc_rdf(job):
    # Create rdf data from the production run
    return

def modify_submit_lammps(filename, statepoint,cores):
    # Modify Submit Scripts
    with open("submit.pbs","r") as f:
        lines = f.readlines()
        lines[1] = "#PBS -N {}{}\n".format(filename, statepoint)
        lines[11] = "mpirun -np {} lmp < {}\n".format(cores, filename)
   with open("submit.pbs","w") as f:
       f.write(lines)
   return

if __name__ == "__main__":
    pr = Project()
    pr.main()
    breakpoint()
