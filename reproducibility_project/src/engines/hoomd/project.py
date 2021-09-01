"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib

import flow
from flow import FlowProject
from flow.environment import DefaultSlurmEnvironment


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
@Project.pre(lambda j: j.doc.get("npt_eq"))
@Project.post(lambda j: j.doc.get("nvt_finished"))
def run_nvt(job):
    """Run a simulation with HOOMD-blue."""
    run_hoomd(job, "nvt", restart=job.isfile("nvt-trajectory.gsd"))


@Project.operation.with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: j.sp.engine == "hoomd")
@Project.post(lambda j: j.doc.get("npt_finished"))
def run_npt(job):
    """Run a simulation with HOOMD-blue."""
    run_hoomd(job, "npt", restart=job.isfile("npt-trajectory.gsd"))


@Project.operation.with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: j.sp.engine == "hoomd")
@Project.pre(lambda j: j.doc.get("npt_finished"))
@Project.post(lambda j: j.doc.get("npt_eq"))
def check_equilibration_npt(job):
    """Run a simulation with HOOMD-blue."""
    job.doc.npt_finished = check_equilibration(job, "npt", "volume")


@Project.operation.with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: j.sp.engine == "hoomd")
@Project.pre(lambda j: j.doc.get("nvt_finished"))
@Project.post(lambda j: j.doc.get("nvt_eq"))
def check_equilibration_nvt(job):
    """Run a simulation with HOOMD-blue."""
    job.doc.nvt_finished = check_equilibration(job, "nvt", "potential_energy")


def run_hoomd(job, method, restart=False):
    """Run a simulation with HOOMD-blue."""
    import foyer
    import gsd.hoomd
    import hoomd
    import hoomd.md
    import numpy as np
    import unyt as u
    from mbuild.formats.gsdwriter import write_gsd
    from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield

    from reproducibility_project.src.molecules.system_builder import (
        construct_system,
    )
    from reproducibility_project.src.utils.forcefields import load_ff

    if method not in ["npt", "nvt"]:
        raise ValueError("Method must be 'nvt' or 'npt'.")

    # This structure will only be used for the initial npt run,
    # but we need forcefield info for all sims.
    # Ignore the vapor box
    # Initializing at high density causes issues, so instead we initialize
    # with box expanded by factor
    filled_box, _ = construct_system(job.sp, scale_liq_box=2)

    ff = load_ff(job.sp.forcefield_name)
    structure = ff.apply(filled_box)

    # ref_distance: 10 angstrom -> 1 nm
    # ref_energy: 1/4.184 kcal/mol -> 1 kJ/mol
    # ref_mass: 0.9999938574 dalton -> 1 amu
    d = 10
    e = 1 / 4.184
    m = 0.9999938574

    snapshot, forcefield, ref_vals = create_hoomd_forcefield(
        structure, ref_distance=d, ref_energy=e, ref_mass=m, r_cut=job.sp.r_cut
    )
    if method == "npt":
        print("NPT")
        if restart:
            print("restart")
            with gsd.hoomd.open(job.fn("trajectory-npt.gsd")) as t:
                snapshot = t[-1]
        else:
            write_gsd(
                structure,
                job.fn("init.gsd"),
                ref_distance=d,
                ref_energy=e,
                ref_mass=m,
            )

    else:
        print("NVT")
        if restart:
            print("restart")
            with gsd.hoomd.open(job.fn("trajectory-nvt.gsd")) as t:
                snapshot = t[-1]
        else:
            # nvt overwrites snapshot information with snapshot from npt run
            with gsd.hoomd.open(job.fn("trajectory-npt.gsd")) as t:
                snapshot = t[-1]

    if restart:
        writemode = "a"
    else:
        writemode = "w"

    device = hoomd.device.auto_select()
    sim = hoomd.Simulation(device=device, seed=job.sp.replica)
    sim.create_state_from_snapshot(snapshot)
    gsd_writer = hoomd.write.GSD(
        filename=job.fn(f"trajectory-{method}.gsd"),
        trigger=hoomd.trigger.Periodic(10000),
        mode=f"{writemode}b",
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
    integrator = hoomd.md.Integrator(dt=dt)
    integrator.forces = forcefield
    # convert temp in K to kJ/mol
    kT = (job.sp.temperature * u.K).to_equivalent("kJ/mol", "thermal").value
    if method == "npt":
        # convert pressure to unit system
        pressure = (job.sp.pressure * u.kPa).to("kJ/(mol*nm**3)").value
        integrator_method = hoomd.md.methods.NPT(
            filter=hoomd.filter.All(),
            kT=kT,
            tau=1000 * dt,
            S=pressure,
            tauS=5000 * dt,
            couple="xyz",
        )
    else:
        integrator_method = hoomd.md.methods.NVT(
            filter=hoomd.filter.All(), kT=kT, tau=1000 * dt
        )
    integrator.methods = [integrator_method]
    sim.operations.integrator = integrator
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT)

    if method == "npt":
        # only run with high tauS if we are starting from scratch
        if not restart:
            sim.run(1e6)
        integrator.tauS = 500 * dt
        print(f"tauS: {sim.operations.integrator.tauS}")
    else:
        if not restart:
            # Shrink step follows this example
            # https://hoomd-blue.readthedocs.io/en/latest/tutorial/
            # 01-Introducing-Molecular-Dynamics/03-Compressing-the-System.html
            ramp = hoomd.variant.Ramp(
                A=0, B=1, t_start=sim.timestep, t_ramp=int(2e4)
            )
            # the target volume should be the average volume from NPT
            target_volume = job.doc.avg_volume
            initial_box = sim.state.box
            L = target_volume ** (1 / 3)
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
    job.doc[f"{method}_finished"] = True


def check_equilibration(job, method, eq_property):
    """Check whether a simulation is equilibrated."""
    import numpy as np

    import reproducibility_project.src.analysis.equilibration as eq

    data = np.genfromtxt(job.fn(f"log-{method}.txt"), names=True)
    prop_data = data[eq_property]
    iseq, _, _, _ = eq.is_equilibrated(prop_data)
    if iseq:
        uncorr, i, g, N = eq.trim_non_equilibrated(prop_data)
        job.doc[f"avg_{eq_property}"] = np.average(uncorr)
        job.doc[f"std_{eq_property}"] = np.std(uncorr)
    job.doc[f"{method}_eq"] = iseq
    return iseq


if __name__ == "__main__":
    pr = Project()
    pr.main()
