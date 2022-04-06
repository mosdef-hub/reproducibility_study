"""Setup for signac, signac-flow, signac-dashboard for this study."""
import datetime
import os
import pathlib

import flow
import numpy as np
from flow import FlowProject
from flow.environment import DefaultSlurmEnvironment

rigid_molecules = ["waterSPCE", "benzeneUA"]


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


# The MOSDEF_PYTHON environment variable is set by running
# echo "export MOSDEF_PYTHON=$(which python)" >> ~/.bashrc
# with the mosdef-study38 conda env active
@Project.operation.with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: j.sp.engine == "hoomd")
@Project.post(lambda j: j.doc.get("shrink_finished"))
def run_shrink(job):
    """Initialize volume for simulation with HOOMD-blue."""
    run_hoomd(job, "shrink")


@Project.operation.with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: j.sp.engine == "hoomd")
@Project.pre(lambda j: j.doc.get("npt_eq"))
@Project.post(lambda j: j.doc.get("nvt_finished"))
def run_nvt(job):
    """Run an NVT simulation with HOOMD-blue."""
    run_hoomd(job, "nvt", restart=job.isfile("trajectory-nvt.gsd"))


@Project.operation.with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: j.sp.engine == "hoomd")
@Project.pre(lambda j: j.doc.get("shrink_finished"))
@Project.post(lambda j: j.doc.get("npt_finished"))
def run_npt(job):
    """Run an NPT simulation with HOOMD-blue."""
    run_hoomd(job, "npt", restart=job.isfile("trajectory-npt.gsd"))


@Project.operation.with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: j.sp.engine == "hoomd")
@Project.pre(lambda j: j.doc.get("npt_finished"))
@Project.post(lambda j: j.doc.get("npt_eq"))
def check_equilibration_npt(job):
    """Check the equilibration of the NPT simulation."""
    job.doc.npt_finished = check_equilibration(job, "npt", "volume")


@Project.operation.with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: j.sp.engine == "hoomd")
@Project.pre(lambda j: j.doc.get("nvt_finished"))
@Project.post(lambda j: j.doc.get("nvt_eq"))
def check_equilibration_nvt(job):
    """Check the equilibration of the NVT simulation."""
    job.doc.nvt_finished = check_equilibration(job, "nvt", "potential_energy")


@Project.operation.with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: j.sp.engine == "hoomd")
@Project.pre(lambda j: j.doc.get("nvt_eq"))
@Project.post(lambda j: j.doc.get("post_processed"))
def post_process(job):
    """Run post-processing on the log files."""
    from shutil import copy

    import numpy.lib.recfunctions as rf
    import unyt as u

    for logfile in [job.fn("log-npt.txt"), job.fn("log-nvt.txt")]:
        # Make a copy, just in case
        backup = f"{logfile}.bkup"
        if not os.path.exists(backup):
            copy(logfile, backup)

        data = np.genfromtxt(logfile, names=True)
        data = clean_data(data)

        system_mass = job.sp.mass * u.amu * job.sp.N_liquid
        volume = data["volume"] * u.nm**3
        density = (system_mass / volume).to("g/cm**3")
        kB = 0.00831446262  # kJ/(mol K)
        pressure_factor = float((1 * u.kJ / u.mol / u.nm ** 3).to("kPa"))

        data = rf.drop_fields(data, ["time_remaining"])
        data = rf.rename_fields(data, {"kinetic_temperature": "temperature"})
        data["temperature"] /= kB
        data["pressure"] *= pressure_factor
        data = rf.append_fields(
            data, "density", np.array(density), usemask=False
        )
        np.savetxt(logfile, data, header=" ".join(data.dtype.names))
    job.doc.post_processed = True


def run_hoomd(job, method, restart=False):
    """Run a simulation with HOOMD-blue."""
    import foyer
    import gsd.hoomd
    import hoomd
    import hoomd.md
    import numpy as np
    import unyt as u
    from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield

    from reproducibility_project.src.molecules.system_builder import (
        construct_system,
        get_molecule,
    )
    from reproducibility_project.src.utils.forcefields import load_ff
    from reproducibility_project.src.utils.rigid import moit

    print(job.sp.molecule)
    if method not in ["npt", "nvt", "shrink"]:
        raise ValueError("Method must be 'nvt', 'npt' or 'shrink'.")

    # For rigid molecules, we need to create an initial snapshot with only the
    # rigid body centers
    # Only the number matters at this point--all other attributes of the
    # snapshot will be adjusted later.
    if job.sp["molecule"] in rigid_molecules:
        isrigid = True
        init_snap = hoomd.Snapshot()
        init_snap.particles.types = ["R"]
        N_mols = job.sp["N_liquid"]
        init_snap.particles.N = N_mols
    else:
        isrigid = False
        init_snap = None

    # This structure will only be used for the initial npt run,
    # but we need forcefield info for all sims.
    # Ignore the vapor box
    # Initializing at high density causes issues, so instead we initialize
    # with box expanded by factor
    filled_box, _ = construct_system(
        job.sp, scale_liq_box=2, fix_orientation=isrigid
    )

    ff = load_ff(job.sp.forcefield_name)
    structure = ff.apply(filled_box)

    # ref_distance: 10 angstrom -> 1 nm
    # ref_energy: 1/4.184 kcal/mol -> 1 kJ/mol
    # ref_mass: 0.9999938574 dalton -> 1 amu
    d = 10
    e = 1 / 4.184
    m = 0.9999938574

    snapshot, forcefield, ref_vals = create_hoomd_forcefield(
        structure,
        ref_distance=d,
        ref_energy=e,
        ref_mass=m,
        r_cut=job.sp.r_cut,
        init_snap=init_snap,
        pppm_kwargs={"Nx": 64, "Ny": 64, "Nz": 64, "order": 7},
    )

    # update the neighborlist exclusions for pentane and benzene
    # these wont be set automatically because their scaling is 0
    # forcefield[0] is LJ pair force and all nlist objects are connected
    if job.sp.molecule == "benzeneUA" or job.sp.molecule == "pentaneUA":
        forcefield[0].nlist.exclusions = ["bond", "1-3", "1-4"]
    if job.sp.molecule == "methaneUA":
        forcefield[0].nlist.exclusions = []

    # Adjust the snapshot rigid bodies
    if isrigid:
        # number of particles per molecule
        N_p = get_molecule(job.sp).n_particles
        mol_inds = [
            np.arange(N_mols + i * N_p, N_mols + i * N_p + N_p)
            for i in range(N_mols)
        ]
        for i, inds in enumerate(mol_inds):
            total_mass = np.sum(snapshot.particles.mass[inds])
            # set the rigid body position at the center of mass
            com = (
                np.sum(
                    snapshot.particles.position[inds]
                    * snapshot.particles.mass[inds, np.newaxis],
                    axis=0,
                )
                / total_mass
            )
            snapshot.particles.position[i] = com
            # set the body attribute for the rigid center and its constituents
            snapshot.particles.body[i] = i
            snapshot.particles.body[inds] = i * np.ones_like(inds)
            # set the rigid body center's mass
            snapshot.particles.mass[i] = np.sum(snapshot.particles.mass[inds])
            # set moment of inertia
            snapshot.particles.moment_inertia[i] = moit(
                snapshot.particles.position[inds],
                snapshot.particles.mass[inds],
                center=com,
            )

        # delete the harmonic bond and angle potentials
        remove = [
            f
            for f in forcefield
            if isinstance(f, hoomd.md.bond.Harmonic)
            or isinstance(f, hoomd.md.angle.Harmonic)
            or isinstance(f, hoomd.md.dihedral.Harmonic)
            or isinstance(f, hoomd.md.dihedral.OPLS)
        ]
        for f in remove:
            forcefield.remove(f)

        # update the neighborlist exclusions for rigid
        forcefield[0].nlist.exclusions = ["body"]

    if job.sp.get("long_range_correction") == "energy_pressure":
        for force in forcefield:
            if isinstance(force, hoomd.md.pair.LJ):
                force.tail_correction = True
                print(f"{force} tail_correction set to {force.tail_correction}")

    device = hoomd.device.auto_select()
    print(f"Running HOOMD version {hoomd.version.version}", flush=True)
    if isinstance(device, hoomd.device.GPU):
        print("HOOMD is running on GPU", flush=True)
        print(f"GPU api version {hoomd.version.gpu_api_version}", flush=True)
    else:
        print("HOOMD is running on CPU", flush=True)

    if method == "shrink":
        print("Starting shrink", flush=True)
        sim = hoomd.Simulation(device=device, seed=job.sp.replica)
        sim.create_state_from_snapshot(snapshot)
        hoomd.write.GSD.write(
            state=sim.state, filename=job.fn("init.gsd"), mode="wb"
        )
        filled_box.save(job.fn("starting_compound.json"))
        initgsd = job.fn("init.gsd")

    elif method == "npt":
        print("Starting NPT", flush=True)
        if restart:
            print("Restarting from last frame of existing gsd", flush=True)
            initgsd = job.fn("trajectory-npt.gsd")
        else:
            # npt overwrites snapshot information with snapshot from shrink run
            initgsd = job.fn("trajectory-shrink.gsd")

    else:
        print("Starting NVT", flush=True)
        if restart:
            print("Restarting from last frame of existing gsd", flush=True)
            initgsd = job.fn("trajectory-nvt.gsd")
        else:
            # nvt overwrites snapshot information with snapshot from npt run
            initgsd = job.fn("trajectory-npt.gsd")

    if restart:
        writemode = "a"
    else:
        writemode = "w"

    sim = hoomd.Simulation(device=device, seed=job.sp.replica)
    sim.create_state_from_gsd(initgsd)

    if isrigid:
        # Because we use the fix_orientation flag with fill box, we can safely
        # assume that all the rigid body constituent particles have the same
        # orientation around the rigid body center. Therefore we can define all
        # rigid bodies using just the first one
        rigid = hoomd.md.constrain.Rigid()
        inds = mol_inds[0]

        r_pos = snapshot.particles.position[0]
        c_pos = snapshot.particles.position[inds]
        c_pos -= r_pos
        c_pos = [tuple(i) for i in c_pos]

        c_types = [
            snapshot.particles.types[i] for i in snapshot.particles.typeid[inds]
        ]

        c_orient = [tuple(i) for i in snapshot.particles.orientation[inds]]

        c_charge = [i for i in snapshot.particles.charge[inds]]

        c_diam = [i for i in snapshot.particles.diameter[inds]]

        rigid.body["R"] = {
            "constituent_types": c_types,
            "positions": c_pos,
            "orientations": c_orient,
            "charges": c_charge,
            "diameters": c_diam,
        }

        for force in forcefield:
            if isinstance(force, hoomd.md.pair.LJ):
                for t in snapshot.particles.types:
                    force.params[("R", t)] = dict(epsilon=0, sigma=0)
                    force.r_cut[("R", t)] = 0

        _all = hoomd.filter.Rigid(("center", "free"))
    else:
        _all = hoomd.filter.All()

    gsd_writer = hoomd.write.GSD(
        filename=job.fn(f"trajectory-{method}.gsd"),
        trigger=hoomd.trigger.Periodic(10000),
        mode=f"{writemode}b",
        dynamic=["momentum"],
    )
    sim.operations.writers.append(gsd_writer)

    logger = hoomd.logging.Logger(categories=["scalar", "string"])
    logger.add(sim, quantities=["timestep", "tps"])
    thermo_props = hoomd.md.compute.ThermodynamicQuantities(filter=_all)
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

    status = Status(sim)
    logger[("Status", "time_remaining")] = (status, "time_remaining", "string")
    table_file = hoomd.write.Table(
        output=open(
            job.fn(f"log-{method}.txt"), mode=f"{writemode}", newline="\n"
        ),
        trigger=hoomd.trigger.Periodic(period=1000),
        logger=logger,
        max_header_len=7,
    )
    sim.operations.writers.append(table_file)

    dt = 0.001
    if isrigid:
        integrator = hoomd.md.Integrator(dt=dt, integrate_rotational_dof=True)
        integrator.rigid = rigid
    else:
        integrator = hoomd.md.Integrator(dt=dt)
    integrator.forces = forcefield

    # convert temp in K to kJ/mol
    kT = (job.sp.temperature * u.K).to_equivalent("kJ/mol", "thermal").value

    # start with high tau and tauS
    tau = 1000 * dt
    tauS = 5000 * dt

    if method == "npt":
        # convert pressure to unit system
        pressure = (job.sp.pressure * u.kPa).to("kJ/(mol*nm**3)").value
        integrator_method = hoomd.md.methods.NPT(
            filter=_all,
            kT=kT,
            tau=tau,
            S=pressure,
            tauS=tauS,
            couple="xyz",
        )
    else:
        # shrink and nvt both use nvt
        integrator_method = hoomd.md.methods.NVT(filter=_all, kT=kT, tau=tau)

    integrator.methods = [integrator_method]
    sim.operations.integrator = integrator
    if not restart:
        sim.state.thermalize_particle_momenta(filter=_all, kT=kT)

    if method == "npt":
        # only run with high tauS if we are starting from scratch
        if not restart:
            steps = 1e6
            print(
                f"""Running {steps} with tau: {integrator.methods[0].tau} and
                tauS: {integrator.methods[0].tauS}"""
            )
            sim.run(steps)
            print("Done")
        integrator.methods[0].tauS = 1000 * dt
        integrator.methods[0].tau = 100 * dt
    else:
        # Shrink and NVT both use NVT method
        if method == "shrink":
            # shrink to the desired box length
            L = job.sp.box_L_liq
            shrink_steps = 1e5
        else:
            # The target volume should be the average volume from NPT
            target_volume = job.doc.avg_volume
            L = target_volume ** (1 / 3)
            shrink_steps = 2e4

        if not restart:
            # Shrink and nvt methods both use shrink step
            # Shrink step follows this example
            # https://hoomd-blue.readthedocs.io/en/latest/tutorial/
            # 01-Introducing-Molecular-Dynamics/03-Compressing-the-System.html
            ramp = hoomd.variant.Ramp(
                A=0, B=1, t_start=sim.timestep, t_ramp=int(shrink_steps)
            )
            initial_box = sim.state.box
            final_box = hoomd.Box(Lx=L, Ly=L, Lz=L)
            box_resize_trigger = hoomd.trigger.Periodic(10)
            box_resize = hoomd.update.BoxResize(
                box1=initial_box,
                box2=final_box,
                variant=ramp,
                trigger=box_resize_trigger,
            )
            sim.operations.updaters.append(box_resize)
            print(
                f"Running shrink {shrink_steps} with tau: {integrator.methods[0].tau}"
            )
            sim.run(shrink_steps + 1)
            print("Done")
            assert sim.state.box == final_box
            sim.operations.updaters.remove(box_resize)
            integrator.methods[0].tau = 100 * dt

    if method != "shrink":
        steps = 5e6
        if method == "npt":
            print(
                f"""Running {steps} with tau: {integrator.methods[0].tau} and
                tauS: {integrator.methods[0].tauS}"""
            )
        else:
            print(f"Running {steps} with tau: {integrator.methods[0].tau}")
        sim.run(steps)
    job.doc[f"{method}_finished"] = True
    print("Finished", flush=True)


def check_equilibration(job, method, eq_property, min_t0=100):
    """Check whether a simulation is equilibrated."""
    import numpy as np
    from pymbar.timeseries import subsampleCorrelatedData as subsample

    import reproducibility_project.src.analysis.equilibration as eq

    data = np.genfromtxt(job.fn(f"log-{method}.txt"), names=True)
    data = clean_data(data)
    prop_data = data[eq_property]
    iseq, _, _, _ = eq.is_equilibrated(prop_data)
    if iseq:
        uncorr, t0, g, N = eq.trim_non_equilibrated(prop_data)
        # Sometimes the trim_non_equilibrated function does not cut off enough
        # of the early fluctuating data
        if t0 < min_t0:
            uncorr = prop_data[min_t0:]
        indices = subsample(uncorr, g=g, conservative=True)
        job.doc[f"avg_{eq_property}"] = np.average(prop_data[indices])
        job.doc[f"std_{eq_property}"] = np.std(prop_data[indices])
    job.doc[f"{method}_eq"] = iseq
    return iseq


def clean_data(data):
    """Delete rows in numpy array which contain nan values.

    The HOOMD Table file writer always writes a header, but for restarted jobs,
    this header is in the middle of the file, which creates nan values.
    """
    return np.delete(data, np.where(np.isnan(data["timestep"]))[0], axis=0)


class Status:
    """Monitor the status of a simulation."""

    def __init__(self, sim):
        self.sim = sim

    @property
    def seconds_remaining(self):
        """Compute the seconds remaining based on the current TPS."""
        try:
            return (self.sim.final_timestep - self.sim.timestep) / self.sim.tps
        except ZeroDivisionError:
            return 0

    @property
    def time_remaining(self):
        """Estimate the time remaining."""
        return str(datetime.timedelta(seconds=self.seconds_remaining))


if __name__ == "__main__":
    pr = Project()
    pr.main()
