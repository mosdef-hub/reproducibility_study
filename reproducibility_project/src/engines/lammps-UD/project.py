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


# ____________________________________________________________________________
"""Setting progress label"""


@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
def lammps_created_box(job):
    return job.isfile("box.lammps")


@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
def lammps_copy_files(job):
    return job.isfile("submit.pbs")


@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
def lammps_minimized_equilibrated_nvt(job):
    return job.isfile("minimized.restart_0")


@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
def lammps_equilibrated_npt(job):
    #TODO: modify the following line to properly checking equlibration
    return job.isfile("equilibrated_npt.restart") and True


@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
def lammps_production(job):
    return job.isfile("production.restart")


@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
def lammps_density_data(job):
    return job.isfile("density.dat")


@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
def lammps_created_gsd(job):
    return job.isfile("prod.gsd")


# _____________________________________________________________________
"""Setting up workflow operation"""


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
@Project.post(lammps_created_box)
@flow.with_job
@flow.cmd
def built_lammps(job):
    # Create a lammps datafile for a specified molecule
    from mbuild.formats.lammpsdata import write_lammpsdata
    from project.src.molecules.system_builder import SystemBuilder

    system = SystemBuilder(job)
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
        system,
        "box.lammps",
        atom_style="full",
        unit_style="real",
        mins=system.get_boundingbox().vectors[0],
        maxs=system.get_boundingbox().vectors[1],
        use_rb_torsions=True,
    )
    return


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
@Project.pre(lammps_created_box)
@Project.post(lammps_copy_files)
@flow.with_job
@flow.cmd
def lammps_cp_files(job):
    lmps_submit_path = "../../src/engine_input/lammps/UD_scripts/submit.slurm"
    lmps_run_path = (
        "../../src/engine_input/lammps/input_scripts/in.*"
    )
    msg = f"cp {lmps_submit_path} {lmps_run_path} ./"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
@Project.pre(lammps_copy_files)
@Project.post(lammps_minimized_equilibrated_nvt)
@flow.with_job
@flow.cmd
def lammps_em_nvt(job):
    in_script_name = "in.minimize"
    msg = f"sbatch submit.slurm {in_script_name} {job.sp.seed} {job.sp.temperature} {job.sp.pressure} {job.sp.cutoff}"
    return msg




@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
@Project.pre(lammps_minimized_equilibrated_nvt)
@Project.post(lammps_equilibrated_npt)
@flow.with_job
@flow.cmd
def lammps_equil_npt(job):
    in_script_name = "in.equil"
    msg = f"sbatch submit.slurm {in_script_name} {job.sp.seed} {job.sp.temperature} {job.sp.pressure} {job.sp.cutoff}"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
@Project.pre(lammps_equilibrated_npt)
@Project.post(lammps_production)
@flow.with_job
@flow.cmd
def lammps_prod(job):
    in_script_name = "in.prod"
    msg = f"sbatch submit.slurm {in_script_name} {job.sp.seed} {job.sp.temperature} {job.sp.pressure} {job.sp.cutoff}"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
@Project.pre(lammps_production)
@flow.with_job
@flow.cmd
def lammps_calc_density(job):
    # Create a density datafile from the production run
    return


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
@Project.pre(lammps_production)
@flow.with_job
@flow.cmd
def lammps_calc_rdf(job):
    # Create rdf data from the production run
    import mbuild as mb
    import MDAnalysis as mda

    traj = mda.coordinates.XTC.XTCReader("prod.xtc")
    top = mda.topology.LAMMPSParser.DATAParser("box.lammps")
    u = mda.Universe(top, traj)
    u.trajectory.next(-1)
    parmed_structure = u.convert_to("PARMED")
    mb.formats.gsdwriter.write_gsd(parmed_structure, "prod.gsd")
    # TODO: Use freud rdf PR to create an RDF from the gsd file
    return


def modify_submit_lammps(filename, statepoint, cores):
    # Modify Submit Scripts
    with open("submit.pbs", "r") as f:
        lines = f.readlines()
        lines[1] = "#SBATCH -J {}{}\n".format(filename, statepoint)
    with open("submit.pbs", "w") as f:
        f.write(lines)
    return


if __name__ == "__main__":
    pr = Project()
    pr.main()
    breakpoint()
