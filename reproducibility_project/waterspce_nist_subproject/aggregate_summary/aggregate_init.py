"""Initialize signac statepoints."""
import itertools
import os

import numpy as np
import signac
import unyt as u
from numpy import ModuleDeprecationWarning


def dict_product(dd):
    """Return the product of the key/values of a dictionary."""
    keys = dd.keys()
    for element in itertools.product(*dd.values()):
        yield dict(zip(keys, element))


molecules = [
    "methaneUA",
    "pentaneUA-flexible_bonds",
    "pentaneUA-constrain_bonds",
    "benzeneUA",
    "waterSPCE",
    "ethanolAA",
]
simulation_engines = [
    "cassandra",
    "mcccs",
    "gomc",
    "gromacs",
    "hoomd",
    "lammps-VU",
    "lammps-UD",
]
md_engines = ["gromacs", "hoomd", "lammps-VU", "lammps-UD"]
mc_engines = ["cassandra", "mcccs", "gomc"]
forcefields = {}
r_cuts = {}
cutoff_styles = ["hard"]
long_range_correction = ["energy_pressure"]
for key in molecules:
    if "UA" in key:
        if "benz" not in key:
            forcefields[key] = "trappe-ua"
        else:
            forcefields[key] = "benzene-ua"
        r_cuts[key] = 14.0 * u.angstrom
    elif "SPCE" in key:
        forcefields[key] = "spce"
        r_cuts[key] = 9 * u.angstrom
    else:
        forcefields[key] = "oplsaa"
        r_cuts[key] = 10 * u.angstrom
g_per_cm3 = u.g / (u.cm * u.cm * u.cm)
masses = {
    "methaneUA": [16.04] * u.amu,
    "pentaneUA-flexible_bonds": [72.15] * u.amu,
    "pentaneUA-constrain_bonds": [72.15] * u.amu,
    "benzeneUA": [78.1118] * u.amu,
    "waterSPCE": [18.0153] * u.amu,
    "ethanolAA": [46.0684] * u.amu,
}
init_density_liq = {
    "methaneUA": [0.3752] * g_per_cm3,
    "pentaneUA-flexible_bonds": [0.5390] * g_per_cm3,
    "pentaneUA-constrain_bonds": [0.5390] * g_per_cm3,
    "benzeneUA": [0.692] * g_per_cm3,
    "waterSPCE": [0.998] * g_per_cm3,
    "ethanolAA": [0.7893] * g_per_cm3,
}
init_density_vap = {
    "methaneUA": [0.0117] * g_per_cm3,
    "pentaneUA-flexible_bonds": [0.019] * g_per_cm3,
    "pentaneUA-constrain_bonds": [0.019] * g_per_cm3,
    "benzeneUA": [None],
    "waterSPCE": [None],
    "ethanolAA": [None],
}
temperatures = {
    "methaneUA": [140.0] * u.K,
    "pentaneUA-flexible_bonds": [372.0] * u.K,
    "pentaneUA-constrain_bonds": [372.0] * u.K,
    "benzeneUA": [450.0] * u.K,
    "waterSPCE": [280.0, 300.0, 320.0] * u.K,
    "ethanolAA": [280.0, 300.0, 320.0] * u.K,
}

pressures = {
    "methaneUA": [1318.0] * u.kPa,
    "pentaneUA-flexible_bonds": [1402.0] * u.kPa,
    "pentaneUA-constrain_bonds": [1402.0] * u.kPa,
    "benzeneUA": [2260.0] * u.kPa,
    "waterSPCE": [101.325, 101.325, 101.325] * u.kPa,
    "ethanolAA": [101.325, 101.325, 101.325] * u.kPa,
}

N_liq_molecules = {
    "methaneUA": [900],
    "pentaneUA-flexible_bonds": [300],
    "pentaneUA-constrain_bonds": [300],
    "benzeneUA": [400],
    "waterSPCE": [1100, 1100, 1100],
    "ethanolAA": [500, 500, 500],
}

N_vap_molecules = {
    "methaneUA": [100],
    "pentaneUA-flexible_bonds": [100],
    "pentaneUA-constrain_bonds": [100],
    "benzeneUA": [None],
    "waterSPCE": [None],
    "ethanolAA": [None],
}

liq_box_lengths = {
    "methaneUA": [39.98] * u.angstrom,
    "pentaneUA-flexible_bonds": [40.55] * u.angstrom,
    "pentaneUA-constrain_bonds": [40.55] * u.angstrom,
    "benzeneUA": [42.17] * u.angstrom,
    "waterSPCE": [32.07] * u.angstrom,
    "ethanolAA": [36.46] * u.angstrom,
}

vap_box_lengths = {
    "methaneUA": [61.06] * u.angstrom,
    "pentaneUA-flexible_bonds": [85.75] * u.angstrom,
    "pentaneUA-constrain_bonds": [85.75] * u.angstrom,
    "benzeneUA": [None],
    "waterSPCE": [None],
    "ethanolAA": [None],
}

ensembles = {
    "methaneUA": ["NPT", "GEMC-NVT"],
    "pentaneUA-flexible_bonds": ["NPT", "GEMC-NVT"],
    "pentaneUA-constrain_bonds": ["NPT", "GEMC-NVT"],
    "benzeneUA": ["NPT", None],
    "waterSPCE": ["NPT", None],
    "ethanolAA": ["NPT", None],
}


pr_root = os.path.join(os.getcwd(), "src")
pr = signac.get_project(pr_root)

# ignore statepoints that are not being tested (gemc only for methane, pentane)
# filter the list of dictionaries
total_statepoints = list()
for molecule in molecules:
    for (
        engine,
        ensemble,
        (temp, press),
        n_liq,
        liq_box_L,
        n_vap,
        vap_box_L,
        (init_liq_den, init_vap_den),
        mass,
        lrc,
        cutoff_style,
    ) in itertools.product(
        simulation_engines,
        ensembles[molecule],
        zip(temperatures[molecule], pressures[molecule]),
        N_liq_molecules[molecule],
        liq_box_lengths[molecule],
        N_vap_molecules[molecule],
        vap_box_lengths[molecule],
        zip(init_density_liq[molecule], init_density_vap[molecule]),
        masses[molecule],
        long_range_correction,
        cutoff_styles,
    ):
        statepoint = {
            "molecule": molecule,
            "engine": engine,
            "temperature": np.round(
                temp.to_value("K"),
                decimals=3,
            ).item(),
            "pressure": np.round(press.to_value("kPa"), decimals=3).item(),
            "ensemble": ensemble if ensemble else None,
            "N_liquid": n_liq,
            "N_vap": n_vap if n_vap else None,
            "box_L_liq": np.round(
                liq_box_L.to_value("nm"),
                decimals=3,
            ).item()
            if liq_box_L
            else None,
            "box_L_vap": np.round(
                vap_box_L.to_value("nm"),
                decimals=3,
            ).item()
            if vap_box_L
            else None,
            "init_liq_den": np.round(
                init_liq_den.to_value(g_per_cm3),
                decimals=3,
            ).item(),
            "init_vap_den": np.round(
                init_vap_den.to_value(g_per_cm3),
                decimals=3,
            ).item()
            if init_vap_den
            else None,
            "mass": np.round(
                mass.to_value("amu"),
                decimals=3,
            ).item(),
            "forcefield_name": forcefields[molecule],
            "cutoff_style": cutoff_style,
            "long_range_correction": lrc,
            "r_cut": np.round(
                r_cuts[molecule].to_value("nm"),
                decimals=3,
            ).item(),
        }
        total_statepoints.append(statepoint)

# print(len(total_statepoints))
indices_to_remove = set()
for i, sp in enumerate(total_statepoints):
    # filter gemc ensembles from md engines
    if sp["ensemble"] == "GEMC-NVT" and sp["engine"] in md_engines:
        indices_to_remove.add(i)

    if sp["ensemble"] == "NPT":
        sp["N_vap"] = None
        sp["box_L_vap"] = None
        sp["init_vap_den"] = None

    if sp["ensemble"] is None:
        indices_to_remove.add(i)

    if (
        sp["engine"] in mc_engines
        and sp["molecule"] == "pentaneUA-flexible_bonds"
    ):
        indices_to_remove.add(i)
    if (
        "lammps" in sp["engine"]
        and sp["molecule"] == "pentaneUA-constrain_bonds"
    ):
        indices_to_remove.add(i)


# now reverse sort the set and remove from inital list
# must be reverse sorted to remove indices on the list in place
# otherwise the list will change size and the indices would change
# print(len(indices_to_remove))
sorted_indicies_to_delete = sorted(list(indices_to_remove), reverse=True)
for idx in sorted_indicies_to_delete:
    del total_statepoints[idx]

for sp in total_statepoints:
    pr.open_job(
        statepoint=sp,
    ).init()
