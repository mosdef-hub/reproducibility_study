"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib

import flow
import numpy as np
from flow import FlowProject
from flow.environment import DefaultSlurmEnvironment

rigid_molecules = ["waterSPCE"]


class Project(FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


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


@Project.label
@Project.pre(lambda j: j.sp.engine == "hoomd")
def OutputThermoData(job):
    """Check if the engine loaded the input files and wrote out thermo data."""
    return job.isfile("log-spe-raw.txt")


@Project.label
@Project.pre(lambda j: j.sp.engine == "hoomd")
def FinishedSPECalc(job):
    """Check if the log-spe.txt has been created."""
    return job.isfile("log-spe.txt")


# The MOSDEF_PYTHON environment variable is set by running
# echo "export MOSDEF_PYTHON=$(which python)" >> ~/.bashrc
# with the mosdef-study38 conda env active
@Project.operation.with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: j.sp.engine == "hoomd")
@Project.post(OutputThermoData)
@flow.with_job
def run_singleframe(job):
    """Create and run initial configurations of the system statepoint."""
    import hoomd
    import hoomd.md
    import mbuild as mb
    import numpy as np
    import unyt as u
    from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield

    from reproducibility_project.src.molecules.system_builder import (
        get_molecule,
    )
    from reproducibility_project.src.utils.forcefields import load_ff
    from reproducibility_project.src.utils.rigid import moit

    molecule = job.sp.molecule
    pr = Project()
    snapshot_directory = (
        pathlib.Path(pr.root_directory()) / "src" / "system_snapshots"
    )
    molecule_filename = molecule + ".json"
    box = mb.load(str(snapshot_directory / molecule_filename))

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

    ff = load_ff(job.sp.forcefield_name)
    structure = ff.apply(box)

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
    print("Snapshot created")

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
        for f in forcefield:
            f.nlist.exclusions = f.nlist.exclusions + ["body"]

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

    sim = hoomd.Simulation(device=device, seed=job.sp.replica)
    sim.create_state_from_snapshot(snapshot)

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

    thermo_props = hoomd.md.compute.ThermodynamicQuantities(filter=_all)
    sim.operations.computes.append(thermo_props)

    dt = 0.001
    if isrigid:
        integrator = hoomd.md.Integrator(dt=dt, integrate_rotational_dof=True)
        integrator.rigid = rigid
    else:
        integrator = hoomd.md.Integrator(dt=dt)
    integrator.forces = forcefield

    # convert temp in K to kJ/mol
    kT = (job.sp.temperature * u.K).to_equivalent("kJ/mol", "thermal").value

    tau = 1000 * dt

    integrator_method = hoomd.md.methods.NVT(filter=_all, kT=kT, tau=tau)

    integrator.methods = [integrator_method]
    sim.operations.integrator = integrator

    sim.run(0)

    labels = []
    values = []
    for force in forcefield:
        labels.append(str(type(force)).split("'")[1])
        values.append(f"{force.energy:.15g}")
        if isinstance(force, hoomd.md.pair.LJ) or isinstance(
            force, hoomd.md.long_range.pppm.Coulomb
        ):
            labels.append(str(type(force)).split("'")[1] + "_tail")
            values.append(f"{force.additional_energy:.15g}")

    for val in ["kinetic_energy", "potential_energy"]:
        labels.append(val)
        values.append(f"{getattr(thermo_props,val):.15g}")

    with open(job.fn("log-spe-raw.txt"), "w") as f:
        f.write(f"{'   '.join(labels)}\n")
        f.write(f"{'   '.join(values)}")

    print("Finished", flush=True)


@Project.operation.with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: j.sp.engine == "hoomd")
@Project.pre(OutputThermoData)
@Project.post(FinishedSPECalc)
@flow.with_job
@flow.cmd
def FormatTextFile(job):
    """Convert simulation engine output to log-spe.txt for data comparisons.

    See README.md for spe_subproject for formatting information.
    """
    import numpy as np
    import unyt as u

    raw_logfile = job.fn("log-spe-raw.txt")
    logfile = job.fn("log-spe.txt")
    data = np.genfromtxt(raw_logfile, names=True)

    headers = [
        "hoomd.md.pair.pair.LJ",
        "hoomd.md.pair.pair.Ewald",
        "hoomd.md.long_range.pppm.Coulomb",
        "hoomd.md.special_pair.LJ",
        "hoomd.md.special_pair.Coulomb",
        "hoomd.md.bond.Harmonic",
        "hoomd.md.angle.Harmonic",
        "hoomd.md.dihedral.OPLS",
        "kinetic_energy",
        "potential_energy",
    ]

    # np.savetxt(logfile,data,header=" ".join(data.dtype.names))


if __name__ == "__main__":
    pr = Project()
    pr.main()
