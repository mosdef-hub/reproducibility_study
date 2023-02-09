"""Summarize spe data and save as csv files."""
import os

import numpy as np
import pandas as pd
import signac
import unyt as u


def K_to_kj_per_mol_conversion(value):
    """Convert energy unit from K to kj/mol."""
    conversion_factor = (1 * u.Kelvin).to_value(
        "kJ/mol", equivalence="thermal"
    )  # conversion_K_to_kj_per_mol
    converted = value * conversion_factor
    return converted


if __name__ == "__main__":
    project = signac.get_project()

    engines = ["lammps-VU", "hoomd", "gromacs", "mcccs", "gomc", "cassandra"]

    molecules = [
        "methaneUA",
        "pentaneUA",
        "benzeneUA",
        "waterSPCE",
        "ethanolAA",
    ]
    properties = [
        "potential_energy",
        "tot_vdw_energy",
        "tail_energy",
        "tot_electrostatics",
        "short_range_electrostatics",
        "long_range_electrostatics",
        "tot_pair_energy",
        "bonds_energy",
        "angles_energy",
        "dihedrals_energy",
        "tot_bonded_energy",
        "intramolecular_energy",
        "intermolecular_energy",
    ]

    spe_data = dict()
    for molecule in molecules:
        spe_data[molecule] = dict()
        for engine in engines:
            spe_data[molecule][engine] = dict()
            for property in properties:
                spe_data[molecule][engine][property] = 0

    # Parse LAMMPS SPE data
    lmp_map = {
        "potential_energy": "potential_energy",
        "vdw_energy": "tot_vdw_energy",
        "coul_energy": "tot_electrostatics",
        "pair_energy": "tot_pair_energy",
        "bonds_energy": "bonds_energy",
        "angles_energy": "angles_energy",
        "dihedrals_energy": "dihedrals_energy",
        "tail_energy": "tail_energy",
        "kspace_energy": "long_range_electrostatics",
    }

    for job in project.find_jobs({"engine": "lammps-VU"}):
        data = np.genfromtxt(f"{job.ws}/log-spe.txt", names=True, delimiter=",")
        for prop in data.dtype.names:
            spe_data[job.sp.molecule][job.sp.engine][lmp_map[prop]] += data[
                prop
            ]

    # Parse HOOMD SPE data
    hoomd_map = {
        "pair_LJ": "tot_pair_energy",
        "pair_LJ_tail": "tail_energy",
        "pair_Ewald": "short_range_electrostatics",
        "pppm_Coulomb": "long_range_electrostatics",
        "special_pair_LJ": "tot_pair_energy",
        "special_pair_Coulomb": "short_range_electrostatics",
        "bond_Harmonic": "bonds_energy",
        "angle_Harmonic": "angles_energy",
        "dihedral_OPLS": "dihedrals_energy",
        "potential_energy": "potential_energy",
    }
    for job in project.find_jobs({"engine": "hoomd"}):
        data = np.genfromtxt(
            f"{job.ws}/log-spe-raw.txt", names=True, delimiter=" "
        )
        for prop in data.dtype.names:
            spe_data[job.sp.molecule][job.sp.engine][hoomd_map[prop]] += data[
                prop
            ]

    # Parse GROMACS SPE data
    for job in project.find_jobs({"engine": "gromacs"}):
        data = np.genfromtxt(
            f"{job.ws}/log-spe-p3m.txt", names=True, delimiter=","
        )
        for prop in data.dtype.names:
            if prop in properties:
                spe_data[job.sp.molecule][job.sp.engine][prop] += data[prop]
            else:
                print(f"Not included {job.sp.engine} {job.sp.molecule} {prop}")

    # Parse MCCCS SPE data
    for job in project.find_jobs({"engine": "mcccs"}):
        data = np.genfromtxt(f"{job.ws}/log-spe.txt", names=True, delimiter=",")
        for prop in data.dtype.names:
            if prop in properties:
                spe_data[job.sp.molecule][job.sp.engine][prop] += data[prop]
            else:
                print(f"Not included {job.sp.engine} {job.sp.molecule} {prop}")

    # Parse Cassandra SPE data
    for job in project.find_jobs({"engine": "cassandra"}):
        data = np.genfromtxt(f"{job.ws}/log-spe.txt", names=True, delimiter=",")
        for prop in data.dtype.names:
            if prop in properties:
                spe_data[job.sp.molecule][job.sp.engine][prop] += data[prop]
            else:
                print(f"Not included {job.sp.engine} {job.sp.molecule} {prop}")

    # Parse GOMC SPE data
    gomc_map = {
        "TOTAL": "potential_energy",
        "INTRA(B)": "intramolecular_energy",
        "INTRA(NB)": "intramolecular_energy",
        "BOND(B)": "bonds_energy",
        "ANGLE(B)": "angles_energy",
        "DIHEDRAL(B)": "dihedrals_energy",
        "INTER(LJ)": "intermolecular_energy",
        "LRC": "tail_energy",
        "TOTAL_ELECT": "total_electrostatics",
    }
    for job in project.find_jobs({"engine": "gomc"}):
        with open(f"{job.ws}/out_production_run.dat") as f:
            data = f.readlines()
        for i, line in enumerate(data):
            if "INITIAL SIMULATION ENERGY" in line:
                for title, value in zip(
                    data[i + 2].split(), data[i + 4].split()
                ):
                    # spe_data[job.id][title] = value
                    if title in gomc_map:
                        spe_data[job.sp.molecule][job.sp.engine][
                            gomc_map[title]
                        ] = K_to_kj_per_mol_conversion(float(value))
                    else:
                        print(f"Not included {job.sp.molecule} {title}")

    # Summary and save to csvs
    if not os.path.isdir("csvs"):
        os.mkdir("csvs")

    summarize_dfs = dict()
    for molecule in spe_data:
        print("Saving", molecule)
        df = pd.DataFrame.from_dict(spe_data[molecule]).transpose()
        df.to_csv(f"csvs/{molecule}_spe.csv")
