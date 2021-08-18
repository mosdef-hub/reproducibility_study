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
    # NOTE: create_hoomd3_forcefield depends on PR
    # https://github.com/mosdef-hub/mbuild/pull/871
    snapshot, forcefield, ref_vals = create_hoomd3_forcefield(
        structure, ref_distance=10, ref_energy=1 / 4.184
    )

    device = hoomd.device.auto_select()
    simulation = hoomd.Simulation(device=device, seed=job.sp.replica)
    simulation.create_state_from_snapshot(snapshot)
    gsd_writer = hoomd.write.GSD(
        filename=job.fn("traj.gsd"),
        trigger=hoomd.trigger.Periodic(10000),
        mode="ab",
    )
    simulation.operations.writers.append(gsd_writer)
    integrator = hoomd.md.Integrator(dt=0.005)
    integrator.forces = forcefield
    nvt = hoomd.md.metehods.NVT(
        filter=hoomd.filter.All(), kT=job.sp.temperature, tau=1.0
    )
    integrator.methods = [nvt]
    simulation.operations.intgrator = integrator
    simulation.state.thermalize_particle_momenta(
        filter=hoomd.filter.All(),
        kT=job.sp.temperature,
    )
    simulation.run(1e6)
    return None


if __name__ == "__main__":
    pr = Project()
    pr.main()
