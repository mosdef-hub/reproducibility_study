# Create a list of mbuild serialized json snapshots to be compared for single point energies (spe)
# Create and return only liquid boxes for all systems
import numpy as np
import unyt as u

from reproducibility_project.src.molecules.system_builder import (
    construct_system,
)


def Serialize_Systems():
    # Create a pseudo statpoint dict to read info to system builder
    molecules = [
        "methaneUA",
        "pentaneUA",
        "benzeneUA",
        "waterSPCE",
        "ethanolAA",
    ]
    N_liq_molecules = {
        "methaneUA": 900,
        "pentaneUA": 300,
        "benzeneUA": 400,
        "waterSPCE": 1100,
        "ethanolAA": 500,
    }
    liq_box_lengths = {
        "methaneUA": [39.98] * u.angstrom,
        "pentaneUA": [40.55] * u.angstrom,
        "benzeneUA": [42.17] * u.angstrom,
        "waterSPCE": [32.07] * u.angstrom,
        "ethanolAA": [36.46] * u.angstrom,
    }
    g_per_cm3 = u.g / (u.cm * u.cm * u.cm)
    init_density_liq = {
        "methaneUA": [0.3752] * g_per_cm3,
        "pentaneUA": [0.5390] * g_per_cm3,
        "benzeneUA": [0.692] * g_per_cm3,
        "waterSPCE": [0.998] * g_per_cm3,
        "ethanolAA": [0.7893] * g_per_cm3,
    }
    masses = {
        "methaneUA": [16.04] * u.amu,
        "pentaneUA": [72.15] * u.amu,
        "benzeneUA": [78.1118] * u.amu,
        "waterSPCE": [18.0153] * u.amu,
        "ethanolAA": [46.0684] * u.amu,
    }

    for molecule in molecules:
        statepoint_info = {
            "molecule": molecule,
            "engine": None,
            "replica": 1,
            "temperature": 300,
            "pressure": 101.325,
            "ensemble": "NPT",
            "N_liquid": N_liq_molecules[molecule],
            "N_vap": None,
            "box_L_liq": np.round(
                liq_box_lengths[molecule].to_value("nm"),
                decimals=3,
            ).item(),
            "box_L_vap": None,
            "init_liq_den": np.round(
                init_density_liq[molecule].to_value(g_per_cm3),
                decimals=3,
            ).item(),
            "init_vap_den": None,
            "mass": np.round(
                masses[molecule].to_value("amu"),
                decimals=3,
            ).item(),
            "forcefield_name": None,
            "cutoff_style": None,
            "long_range_correction": None,
            "r_cut": None,
        }
        print(f"Serializing {molecule} snapshots for SPE calculations.")
        system = construct_system(statepoint_info)[0]
        system.save(molecule + ".json", overwrite=False)
        print("_____________________________________\n\n")


if __name__ == "__main__":
    Serialize_Systems()
