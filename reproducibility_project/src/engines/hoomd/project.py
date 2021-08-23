"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib

import flow
from flow import FlowProject
from flow.environment import DefaultSlurmEnvironment

from reproducibility_project.src.utils.forcefields import load_ff


class Project(FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


class Fry(DefaultSlurmEnvironment):
    """Subclass of DefaultSlurmEnvironment for BSU's Fry cluster."""

    hostname_pattern = "fry.boisestate.edu"
    template = "fry.sh"

    @classmethod
    def add_args(cls, parser):
        """Add command line arguments to the submit call."""
        parser.add_argument(
            "--partition",
            default="batch",
            help="Specify the partition to submit to.",
        )
        parser.add_argument("--nodelist", help="Specify the node to submit to.")


# This environment variable is set by running
# echo "export MOSDEF_PYTHON=$(which python)" >> ~/.bashrc
# with the mosdef-study38 conda env active
@Project.operation.with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: j.sp.engine == "hoomd")
@Project.post(lambda j: j.doc.get("finished"))
def run_hoomd(job):
    """Run a simulation with HOOMD-blue."""
    import foyer
    import hoomd
    import hoomd.md
    import unyt as u
    from mbuild.formats.gsdwriter import write_gsd
    from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield

    from reproducibility_project.src.molecules.system_builder import (
        construct_system,
    )

    # temporary hack until benzene and ethanol are added
    try:
        # Ignore the vapor box
        # Initialize with box expanded by factor of 5
        # We will shrink it later
        filled_box, _ = construct_system(
            job.sp, scale_liq_box=5, scale_vap_box=5
        )
    except AttributeError:
        return

    ff = load_ff(job.sp.forcefield_name)
    structure = ff.apply(filled_box)

    # ref_distance: 10 angstrom -> 1 nm
    # ref_energy: 1/4.184 kcal/mol -> 1 kJ/mol
    # ref_mass: 0.9999938574 dalton -> 1 amu
    d = 10
    e = 1 / 4.184
    m = 0.9999938574
    write_gsd(
        structure, job.fn("init.gsd"), ref_distance=d, ref_energy=e, ref_mass=m
    )

    snapshot, forcefield, ref_vals = create_hoomd_forcefield(
        structure, ref_distance=d, ref_energy=e, ref_mass=m
    )

    device = hoomd.device.auto_select()
    sim = hoomd.Simulation(device=device, seed=job.sp.replica)
    sim.create_state_from_snapshot(snapshot)
    gsd_writer = hoomd.write.GSD(
        filename=job.fn("trajectory.gsd"),
        trigger=hoomd.trigger.Periodic(10000),
        mode="ab",
        dynamic=["momentum"],
    )
    sim.operations.writers.append(gsd_writer)

    logger = hoomd.logging.Logger(categories=["scalar"])
    logger.add(sim, quantities=["timestep", "tps"])
    thermo_props = hoomd.md.compute.ThermodynamicQuantities(
        filter=hoomd.filter.All()
    )
    sim.operations.computes.append(thermo_props)
    logger.add(
        thermo_props,
        quantities=[
            "kinetic_energy",
            "potential_energy",
            "pressure",
            "kinetic_temperature",
            "volume",
        ],
    )
    file = open(job.fn("log.txt"), mode="a", newline="\n")
    table_file = hoomd.write.Table(
        output=file,
        trigger=hoomd.trigger.Periodic(period=5000),
        logger=logger,
        max_header_len=7,
    )
    sim.operations.writers.append(table_file)

    integrator = hoomd.md.Integrator(dt=0.005)
    integrator.forces = forcefield
    # convert temp in K to kJ/mol
    kT = (job.sp.temperature * u.K).to_equivalent("kJ/mol", "thermal").value
    nvt = hoomd.md.methods.NVT(filter=hoomd.filter.All(), kT=kT, tau=1.0)
    integrator.methods = [nvt]
    sim.operations.integrator = integrator
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT)

    # Shrink step follows this example
    # https://hoomd-blue.readthedocs.io/en/latest/tutorial/
    # 01-Introducing-Molecular-Dynamics/03-Compressing-the-System.html
    ramp = hoomd.variant.Ramp(A=0, B=1, t_start=sim.timestep, t_ramp=int(2e4))
    initial_box = sim.state.box
    L = job.sp.box_L_liq
    final_box = hoomd.Box(Lx=L, Ly=L, Lz=L)
    box_resize_trigger = hoomd.trigger.Periodic(10)
    box_resize = hoomd.update.BoxResize(
        box1=initial_box,
        box2=final_box,
        variant=ramp,
        trigger=box_resize_trigger,
    )
    sim.operations.updaters.append(box_resize)
    sim.run(2e4 + 1)
    assert sim.state.box == final_box
    sim.operations.updaters.remove(box_resize)

    sim.run(1e6)
    job.doc.finished = True


if __name__ == "__main__":
    pr = Project()
    pr.main()
