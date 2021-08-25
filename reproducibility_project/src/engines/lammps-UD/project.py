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
    """Check if the lammps simulation box has been created for the job."""
    return job.isfile("box.lammps")


@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
def lammps_copy_files(job):
    """Check if the submission scripts have been copied over for the job."""
    return job.isfile("submit.pbs")


@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
def lammps_minimized_equilibrated_nvt(job):
    """Check if the lammps minimization step has run for the job."""
    return job.isfile("minimized.restart_0")


@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
def lammps_equilibrated_npt(job):
    """Check if the lammps equilibration step has run and passed is_equilibrated for the job."""
    return job.isfile("equilibrated_npt.restart") and True


@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
def lammps_production(job):
    """Check if the lammps production step has run for the job."""
    return job.isfile("production.restart")


@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
def lammps_density_data(job):
    """Check if lammps has output density information for the job."""
    return job.isfile("density.dat")


@Project.label
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
def lammps_created_gsd(job):
    """Check if the mdtraj has converted the production to a gsd trajectory for the job."""
    return job.isfile("prod.gsd")


# _____________________________________________________________________
"""Setting up workflow operation"""


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
@Project.post(lammps_created_box)
@flow.with_job
@flow.cmd
def built_lammps(job):
    """Create initial configurations of the system statepoint."""
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
    """Copy over run files for lammps and the SLURM scheduler."""
    lmps_submit_path = "../../src/engine_input/lammps/UD_scripts/submit.slurm"
    lmps_run_path = "../../src/engine_input/lammps/input_scripts/in.*"
    msg = f"cp {lmps_submit_path} {lmps_run_path} ./"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
@Project.pre(lammps_copy_files)
@Project.post(lammps_minimized_equilibrated_nvt)
@flow.with_job
@flow.cmd
def lammps_em_nvt(job):
    """Run energy minimization and nvt ensemble."""
    in_script_name = "in.minimize"
    modify_submit_lammps(in_script_name, job.sp)
    msg = f"sbatch submit.slurm {in_script_name} {job.sp.replica} {job.sp.temperature} {job.sp.pressure} {job.sp.cutoff}"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
@Project.pre(lammps_minimized_equilibrated_nvt)
@Project.post(lammps_equilibrated_npt)
@flow.with_job
@flow.cmd
def lammps_equil_npt(job):
    """Run npt ensemble equilibration."""
    in_script_name = "in.equil"
    modify_submit_lammps(in_script_name, job.sp)
    msg = f"sbatch submit.slurm {in_script_name} {job.sp.replica} {job.sp.temperature} {job.sp.pressure} {job.sp.cutoff}"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
@Project.pre(lammps_equilibrated_npt)
@Project.post(lammps_production)
@flow.with_job
@flow.cmd
def lammps_prod(job):
    """Run npt ensemble production."""
    in_script_name = "in.prod"
    modify_submit_lammps(in_script_name, job.sp)
    msg = f"sbatch submit.slurm {in_script_name} {job.sp.replica} {job.sp.temperature} {job.sp.pressure} {job.sp.cutoff}"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
@Project.pre(lammps_production)
@flow.with_job
@flow.cmd
def lammps_calc_density(job):
    """Create a density text file."""
    return


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "lammps-UD")
@Project.pre(lammps_production)
@flow.with_job
@flow.cmd
def lammps_calc_rdf(job):
    """Create an rdf from the gsd file using Freud analysis scripts."""
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


# TODO: modify this for your purpose
def modify_submit_lammps(filename, statepoint):
    """Modify the submission scripts to include the job and simulation type in the header."""
    with open("submit.slurm", "r") as f:
        lines = f.readlines()
        lines[1] = "#SBATCH -J {}{}\n".format(filename, statepoint)
    with open("submit.slurm", "w") as f:
        f.write(lines)
    return


if __name__ == "__main__":
    pr = Project()
    pr.main()
    breakpoint()
