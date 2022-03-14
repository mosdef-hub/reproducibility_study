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
    import foyer
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
        # forcefield[0] is LJ pair force and all nlist objects are connected
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

    print(job.sp.molecule)
    print("foyer version: ", foyer.__version__)
    print("nlist exclusions: ", forcefield[0].nlist.exclusions)

    sim.run(0)

    labels = []
    values = []
    for force in forcefield:
        label = "_".join(str(type(force)).split("'")[1].split(".")[-2:])
        labels.append(label)
        values.append(f"{force.energy:.15g}")
        if isinstance(force, hoomd.md.pair.LJ):
            labels.append(label + "_tail")
            values.append(f"{force.additional_energy:.15g}")

    labels.append("potential_energy")
    values.append(f"{thermo_props.potential_energy:.15g}")

    with open(job.fn("log-spe-raw.txt"), "w") as f:
        f.write(f"{' '.join(labels)}\n")
        f.write(f"{' '.join(values)}")

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
        "pair_LJ",
        "pair_LJ_tail",
        "pair_Ewald",
        "pppm_Coulomb",
        "special_pair_LJ",
        "special_pair_Coulomb",
        "bond_Harmonic",
        "angle_Harmonic",
        "dihedral_OPLS",
        "potential_energy",
    ]

    data_dict = {}
    for key in headers:
        try:
            val = data[key]
        except ValueError:
            val = None
        data_dict[key] = float(val)

    lj = data_dict["pair_LJ"] + data_dict["special_pair_LJ"]
    tail = data_dict["pair_LJ_tail"]
    ewald = data_dict["pair_Ewald"]
    coulomb = data_dict["pppm_Coulomb"] + data_dict["special_pair_Coulomb"]

    new_data_dict = {
        "potential_energy": data_dict["potential_energy"],
        "tot_vdw_energy": lj,
        "tail_energy": tail,
        "short_range_electrostatics": ewald,
        "long_range_electrostatics": coulomb,
        "tot_pair_energy": lj + ewald + coulomb,
        "bonds_energy": data_dict["bond_Harmonic"],
        "angles_energy": data_dict["angle_Harmonic"],
        "dihedrals_energy": data_dict["dihedral_OPLS"],
        "tot_bonded_energy": data_dict["bond_Harmonic"]
        + data_dict["angle_Harmonic"]
        + data_dict["dihedral_OPLS"],
        "tot_electrostatics": ewald + coulomb,
        "intramolecular_energy": None,
        "intermolecular_energy": None,
    }

    # dicts are ordered in python3.6+
    with open(logfile, "w") as f:
        f.write(" ".join(new_data_dict.keys()) + "\n")
        f.write(" ".join([f"{i:.15g}" for i in new_data_dict.values()]))

    print("Finished", flush=True)


if __name__ == "__main__":
    pr = Project()
    pr.main()
