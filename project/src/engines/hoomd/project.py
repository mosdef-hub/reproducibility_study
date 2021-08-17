"""Setup for signac, signac-flow, signac-dashboard for this study."""
# import foyer
import os
import pathlib

import flow
import foyer
from flow import environments


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "hoomd")
def run_hoomd(job):
    """Run a simulation with HOOMD-blue."""
    import hoomd
    import hoomd.md
    from mbuild.formats.gsdwriter import write_gsd
    from mbuild.formats.hoomd3_simulation import create_hoomd3_forcefield

    filled_box = get_system(job)
    # ff = foyer.Forcefield(job._project.ff_fn)
    structure = ff.apply(filled_box)

    write_gsd(structure, job.fn("init.gsd"), ref_distance=rd, ref_energy=re)
    # ref_distance: angstrom -> nm
    # ref_energy: kcal/mol -> kJ/mol
    # ref_mass: dalton -> amu
    snapshot, forcefield, ref_vals = create_hoomd3_forcefield(
        structure, ref_distance=0.1, ref_energy=1 / 4.184, ref_mass=0.9999938574
    )

    device = hoomd.device.auto_select()
    sim = hoomd.Simulation(device=device, seed=job.sp.replica)
    sim.create_state_from_snapshot(snapshot)
    gsd_writer = hoomd.write.GSD(
        filename=job.fn("trajectory.gsd"),
        trigger=hoomd.trigger.Periodic(10000),
        mode="ab",
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
    file = open("log.txt", mode="x", newline="\n")
    table_file = hoomd.write.Table(
        output=file,
        trigger=hoomd.trigger.Periodic(period=5000),
        logger=logger,
        max_header_len=7,
    )
    sim.operations.writers.append(table_file)

    integrator = hoomd.md.Integrator(dt=0.005)
    integrator.forces = forcefield
    nvt = hoomd.md.metehods.NVT(
        filter=hoomd.filter.All(), kT=job.sp.temperature, tau=1.0
    )
    integrator.methods = [nvt]
    sim.operations.intgrator = integrator
    sim.state.thermalize_particle_momenta(
        filter=hoomd.filter.All(),
        kT=job.sp.temperature,
    )
    sim.run(1e6)


if __name__ == "__main__":
    pr = Project()
    pr.main()
