"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib
import sys

import flow
import numpy as np
from flow import FlowProject
from flow.environment import DefaultSlurmEnvironment

rigid_molecules = ["waterSPCE"]  # , "benzeneUA"]


class Project(FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


class RahmanHOOMD(DefaultSlurmEnvironment):
    """Subclass of DefaultSlurmEnvironment for VU's Rahman cluster."""

    template = "rahman_hoomd.sh"

    @classmethod
    def add_args(cls, parser):
        """Add command line arguments to the submit call."""
        parser.add_argument(
            "--walltime",
            type=float,
            default=96,
            help="Walltime for this submission",
        )


@Project.label
@Project.pre(lambda j: "hoomd" in j.sp.engine)
def OutputThermoData(job):
    """Check if the engine loaded the input files and wrote out thermo data."""
    return job.isfile("log-spe-raw.txt")


@Project.label
@Project.pre(lambda j: "hoomd" in j.sp.engine)
def FinishedSPECalc(job):
    """Check if the log-spe.txt has been created."""
    return job.isfile("log-spe.txt")


# The MOSDEF_PYTHON environment variable is set by running
# echo "export MOSDEF_PYTHON=$(which python)" >> ~/.bashrc
# with the mosdef-study38 conda env active
@Project.operation  # .with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: "hoomd" in j.sp.engine)
@Project.post(OutputThermoData)
@flow.with_job
def run_singleframe(job):
    """Create and run initial configurations of the system statepoint."""
    import foyer

    # import git
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

    # repo = git.Repo(search_parent_directories=True)
    # sha = repo.head.object.hexsha
    molecule = job.sp.molecule
    print(job.sp.molecule)
    # print(f"git commit: {sha}\n")

    pr = Project()
    snapshot_directory = (
        pathlib.Path(pr.root_directory()) / "src" / "system_snapshots"
    )
    molecule_filename = molecule + ".json"
    box = mb.load(str(snapshot_directory / molecule_filename))
    print(box)
    print(box.name)
    print(list(box.particles())[0])

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
    print(ff)
    print(ff.name)
    print(ff._atomTypes)
    structure = ff.apply(box, use_residue_map=False)

    # ref_distance: 10 angstrom -> 1 nm
    # ref_energy: 1/4.184 kcal/mol -> 1 kJ/mol
    # ref_mass: 0.9999938574 dalton -> 1 amu
    d = 10
    e = 1 / 4.184
    m = 0.9999938574

    if job.sp.molecule == "ethanolAA":
        pppmDict = {"Nx": 24, "Ny": 24, "Nz": 24, "order": 5}
    if job.sp.molecule == "waterSPCE":
        pppmDict = {"Nx": 32, "Ny": 32, "Nz": 32, "order": 5}
    else:
        pppmDict = {"Nx": 64, "Ny": 64, "Nz": 64, "order": 7}

    snapshot, forcefield, ref_vals = create_hoomd_forcefield(
        structure,
        ref_distance=d,
        ref_energy=e,
        ref_mass=m,
        r_cut=job.sp.r_cut,
        init_snap=init_snap,
        pppm_kwargs=pppmDict,
    )
    print("Snapshot created")
    print(f"box: {snapshot.configuration.box}")

    # Adjust the snapshot rigid bodies
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
            or isinstance(f, hoomd.md.dihedral.OPLS)
        ]
        for f in remove:
            forcefield.remove(f)

    if isrigid:
        # Because we use the fix_orientation flag with fill box, we can safely
        # assume that all the rigid body constituent particles have the same
        # orientation around the rigid body center. Therefore we can define all
        # rigid bodies using just the first one
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
            # "charges": c_charge,
            # "diameters": c_diam,
        }

        for force in forcefield:
            if isinstance(force, hoomd.md.pair.LJ):
                for t in snapshot.particles.types:
                    force.params[("R", t)] = dict(epsilon=0, sigma=0)
                    force.r_cut[("R", t)] = 0

    # update the neighborlist exclusions for pentane and methane
    # pentane's wont be set automatically because the scaling is 0
    # and the default (bond, 1-3) is unecessary and raises a warning for methane
    if job.sp.molecule == "pentaneUA":
        forcefield[0].nlist.exclusions = ["bond", "1-3", "1-4"]
    if job.sp.molecule == "benzeneUA":
        forcefield[0].nlist.exclusions = [
            "bond",
            "1-3",
            "1-4",
            "body",
            "angle",
            "dihedral",
        ]
    elif job.sp.molecule == "methaneUA":
        forcefield[0].nlist.exclusions = []

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
    tauS = 1000 * dt
    pressure = (job.sp.pressure * u.kPa).to("kJ/(mol*nm**3)").value
    integrator_method = hoomd.md.methods.ConstantPressure(
        thermostat=hoomd.md.methods.thermostats.Bussi(kT),
        filter=_all,
        S=pressure,
        tauS=tauS,
        couple="xyz",
    )

    integrator.methods = [integrator_method]
    sim.operations.integrator = integrator

    print("mbuild version: ", mb.__version__)
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


@Project.operation  # .with_directives({"executable": "$MOSDEF_PYTHON", "ngpu": 1})
@Project.pre(lambda j: "hoomd" in j.sp.engine)
@Project.pre(OutputThermoData)
@Project.post(FinishedSPECalc)
@flow.with_job
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
            val = 0
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
        "intramolecular_energy": 0,
        "intermolecular_energy": 0,
    }

    # dicts are ordered in python3.6+
    with open(logfile, "w") as f:
        f.write(" ".join(new_data_dict.keys()) + "\n")
        f.write(" ".join([f"{i:.15g}" for i in new_data_dict.values()]))

    print("Finished", flush=True)


if __name__ == "__main__":
    pr = Project()
    pr.main()
