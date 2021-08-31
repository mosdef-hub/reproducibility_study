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
@Project.pre(lambda j: j.doc.get("npt_finished"))
@Project.post(lambda j: j.doc.get("nvt_finished"))
def run_nvt(job):
    """Run a simulation with HOOMD-blue."""
    run_hoomd(job, "nvt")


@Project.operation.with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: j.sp.engine == "hoomd")
@Project.post(lambda j: j.doc.get("npt_finished"))
def run_npt(job):
    """Run a simulation with HOOMD-blue."""
    run_hoomd(job, "npt")


def run_hoomd(job, method):
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

    # This structure will only be used for npt, but we need ff info for both
    # Ignore the vapor box
    # Initializing at high density causes issues, so instead we initialize
    # with box expanded by factor
    filled_box, _ = construct_system(job.sp, scale_liq_box=2.5)

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

        write_gsd(
            structure,
            job.fn("init.gsd"),
            ref_distance=d,
            ref_energy=e,
            ref_mass=m,
        )

    else:
        print("NVT")
        # nvt overwrites snapshot information with snapshot from npt run
        with gsd.hoomd.open(job.fn("trajectory-npt.gsd")) as t:
            snapshot = t[-1]


    device = hoomd.device.auto_select()
    sim = hoomd.Simulation(device=device, seed=job.sp.replica)
    sim.create_state_from_snapshot(snapshot)
    gsd_writer = hoomd.write.GSD(
        filename=job.fn(f"trajectory-{method}.gsd"),
        trigger=hoomd.trigger.Periodic(10000),
        mode="wb",
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
        output=open(job.fn(f"log-{method}.txt"), mode="a", newline="\n")
        trigger=hoomd.trigger.Periodic(period=1000),
        logger=logger,
        max_header_len=7,
    )
    sim.operations.writers.append(table_file)

    integrator = hoomd.md.Integrator(dt=0.005)
    integrator.forces = forcefield
    # convert temp in K to kJ/mol
    kT = (job.sp.temperature * u.K).to_equivalent("kJ/mol", "thermal").value
    if method == "npt":
        # convert pressure to unit system
        pressure = (job.sp.pressure * u.kPa).to("kJ/(mol*nm**3)").value
        integrator_method = hoomd.md.methods.NPT(
            filter=hoomd.filter.All(),
            kT=kT,
            tau=1.0,
            S=pressure,
            tauS=1.0,
            couple="xyz",
        )
    else:
        integrator_method = hoomd.md.methods.NVT(
            filter=hoomd.filter.All(), kT=kT, tau=1.0
        )
    integrator.methods = [integrator_method]
    sim.operations.integrator = integrator
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT)

    sim.run(1e6)
    if method == "npt":
        job.doc.npt_finished = True
    else:
        job.doc.nvt_finished = True


if __name__ == "__main__":
    pr = Project()
    pr.main()
