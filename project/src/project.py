"""Setup for signac, signac-flow, signac-dashboard for this study."""
import flow
from flow import environments
import foyer


class Project(flow.FlowProject):
    pass


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine=='hoomd')
def run_hoomd(job):
    from mbuild.formats.gsdwriter import write_gsd
    import hoomd
    import hoomd.md
    filled_box = get_system(job)z
    ff_file = get_ff_file(job)
    ff = foyer.Forcefield(ff_file)
    structure = ff.apply(filled_box)

    write_gsd(structure, job.fn('init.gsd'), ref_distance=rd, ref_energy=re)
    snapshot, forcefield, ref_vals = create_hoomd3_forcefield(structure,
            ref_distance=10, ref_energy=1/4.184)

    device = hoomd.device.auto_select()
    simulation = hoomd.Simulation(device=device, seed=job.sp.replica)
    simulation.create_state_from_snapshot(snapshot)
    gsd_writer = hoomd.write.GSD(filename=job.fn('traj.gsd'),
            trigger=hoomd.trigger.Periodic(10000),
            mode='ab')
    simulation.operations.writers.append(gsd_writer)
    integrator = hoomd.md.Integrator(dt=0.005)
    integrator.forces = forcefield
    nvt = hoomd.md.metehods.NVT(
            filter=hoomd.filter.All(),
            kT=job.sp.temperature,
            tau=1.0
    )
    integrator.methods = [nvt]
    simulation.operations.intgrator = integrator
    simulation.state.thermalize_particle_momenta(
            filter=hoomd.filter.All(),
            kT=job.sp.temperature,
    )
    simulation.run(1e6)
    return None


if __name__ == '__main__':
    pr = Project()
    pr.main()
