"""Setup for signac, signac-flow, signac-dashboard for this study."""
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


@Project.operation
@Project.pre(lambda j: j.sp.engine == "hoomd")
def run_hoomd(job):
    """Run a simulation with HOOMD-blue."""
    import sys

    import foyer
    import hoomd
    import hoomd.md
    import unyt as u
    from mbuild.formats.gsdwriter import write_gsd
    from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield

    sys.path.append(Project().root_directory())
    from src.molecules.system_builder import construct_system

    # temporary hack until benzene and ethanol are added
    try:
        # Ignore the vapor box
        filled_box, _ = construct_system(job.sp)
    except AttributeError:
        return

    if job.sp.forcefield_name == "trappe-ua":
        ff = foyer.forcefields.load_TRAPPE_UA()
    elif job.sp.forcefield_name == "benzene-ua":
        ff = foyer.forcefields.Forcefield("src/xmls/benzene_trappe-ua_like.xml")
    elif job.sp.forcefield_name == "spce":
        ff = foyer.forcefields.Forcefield("src/xmls/spce.xml")
    elif job.sp.forcefield_name == "oplsaa":
        ff = foyer.forcefields.load_OPLSAA()
    else:
        raise ValueError(f"No forcefield defined for {job.sp.forcefield_name}")

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
    file = open("log.txt", mode="a", newline="\n")
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
    sim.run(1e6)


if __name__ == "__main__":
    pr = Project()
    pr.main()
