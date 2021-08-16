"""Initialize signac statepoints."""
import itertools
import os

import signac
import unyt as u
from numpy import ModuleDeprecationWarning


def dict_product(dd):
    """Return the product of the key/values of a dictionary."""
    keys = dd.keys()
    for element in itertools.product(*dd.values()):
        yield dict(zip(keys, element))


molecules = ["methaneUA", "pentaneUA", "benzeneUA", "waterSPC/E", "ethanolAA"]
replicas = range(5)
simulation_engines = [
    "cassandra",
    "mcccs",
    "gomc",
    "gromacs",
    "hoomd",
    "lammps",
]
md_engines = ["gromacs", "hoomd", "lammps"]
mc_engines = ["cassandra", "mcccs", "gomc"]
forcefields = {}
r_cuts = {}
cutoff_styles = ["hard"]
for key in molecules:
    if "UA" in key:
        if "benz" not in key:
            forcefields[key] = "Trappe_UA"
        else:
            forcefields[key] = "benzene_ua"
        r_cuts[key] = 14.0 * u.angstrom
    elif "SPC/E" in key:
        forcefields[key] = "spce"
        r_cuts[key] = 9 * u.angstrom
    else:
        forcefields[key] = "oplsaa"
        r_cuts[key] = 10 * u.angstrom
g_per_mol = u.g / u.mol
g_per_cm3 = u.g / (u.cm * u.cm * u.cm)
masses = {
    "methaneUA": [16.04] * g_per_mol,
    "pentaneUA": [72.15] * g_per_mol,
    "benzeneUA": [78.1118] * g_per_mol,
    "waterSPC/E": [18.0153] * g_per_mol,
    "ethanolAA": [46.0684] * g_per_mol,
}
init_density_liq = {
    "methaneUA": [0.3752] * g_per_cm3,
    "pentaneUA": [0.5390] * g_per_cm3,
    "benzeneUA": [0.692] * g_per_cm3,
    "waterSPC/E": [0.998] * g_per_cm3,
    "ethanolAA": [0.7893] * g_per_cm3,
}
init_density_vap = {
    "methaneUA": [0.0117] * g_per_cm3,
    "pentaneUA": [0.019] * g_per_cm3,
    "benzeneUA": [None],
    "waterSPC/E": [None],
    "ethanolAA": [None],
}
temperatures = {
    "methaneUA": [140.0] * u.K,
    "pentaneUA": [372.0] * u.K,
    "benzeneUA": [450.0] * u.K,
    "waterSPC/E": [280.0, 300.0, 320.0] * u.K,
    "ethanolAA": [280.0, 300.0, 320.0] * u.K,
}

pressures = {
    "methaneUA": [1318.0] * u.kPa,
    "pentaneUA": [1402.0] * u.kPa,
    "benzeneUA": [2260.0] * u.kPa,
    "waterSPC/E": [101.325, 101.325, 101.325] * u.kPa,
    "ethanolAA": [101.325, 101.325, 101.325] * u.kPa,
}

N_liq_molecules = {
    "methaneUA": [900],
    "pentaneUA": [300],
    "benzeneUA": [400],
    "waterSPC/E": [2000, 2000, 2000],
    "ethanolAA": [700, 700, 700],
}

N_vap_molecules = {
    "methaneUA": [100],
    "pentaneUA": [100],
    "benzeneUA": [None],
    "waterSPC/E": [None],
    "ethanolAA": [None],
}

liq_box_lengths = {
    "methaneUA": [39.98] * u.angstrom,
    "pentaneUA": [40.55] * u.angstrom,
    "benzeneUA": [42.17] * u.angstrom,
    "waterSPC/E": [39.14] * u.angstrom,
    "ethanolAA": [40.79] * u.angstrom,
}

vap_box_lengths = {
    "methaneUA": [61.06] * u.angstrom,
    "pentaneUA": [85.75] * u.angstrom,
    "benzeneUA": [None],
    "waterSPC/E": [None],
    "ethanolAA": [None],
}

ensembles = {
    "methaneUA": ["NPT", "GEMC-NVT"],
    "pentaneUA": ["NPT", "GEMC-NVT"],
    "benzeneUA": ["NPT", None],
    "waterSPC/E": ["NPT", None],
    "ethanolAA": ["NPT", None],
}


pr_root = os.path.join(os.getcwd(), "src")
pr = signac.get_project(pr_root)

# ignore statepoints that are not being tested (gemc only for methane, pentane)
# filter the list of dictionaries
total_statepoints = list()
for molecule in molecules:
    for engine in simulation_engines:
        for ensemble in ensembles[molecule]:
            for temp, press in zip(temperatures[molecule], pressures[molecule]):
                for n_liq in N_liq_molecules[molecule]:
                    for liq_box_L in liq_box_lengths[molecule]:
                        for n_vap in N_vap_molecules[molecule]:
                            for vap_box_L in vap_box_lengths[molecule]:
                                for init_liq_den, init_vap_den in zip(
                                    init_density_liq[molecule],
                                    init_density_vap[molecule],
                                ):
                                    for mass in masses[molecule]:
                                        for cutoff_style in cutoff_styles:
                                            for replica in replicas:
                                                statepoint = {
                                                    "molecule": molecule,
                                                    "engine": engine,
                                                    "replica": replica,
                                                    "temperature": temp.to_value(
                                                        "K"
                                                    ),
                                                    "pressure": press.to_value(
                                                        "kPa"
                                                    ),
                                                    "ensemble": ensemble
                                                    if ensemble
                                                    else None,
                                                    "N_liquid": n_liq,
                                                    "N_vap": n_vap
                                                    if n_vap
                                                    else None,
                                                    "box_L_liq": liq_box_L.to_value(
                                                        "nm"
                                                    )
                                                    if liq_box_L
                                                    else None,
                                                    "box_L_vap": vap_box_L.to_value(
                                                        "nm"
                                                    )
                                                    if vap_box_L
                                                    else None,
                                                    "init_liq_den": init_liq_den.to_value(
                                                        g_per_cm3
                                                    ),
                                                    "init_vap_den": init_vap_den.to_value(
                                                        g_per_cm3
                                                    )
                                                    if init_vap_den
                                                    else None,
                                                    "mass": mass.to_value(
                                                        "g/mol"
                                                    ),
                                                    "forcefield_name": forcefields[
                                                        molecule
                                                    ],
                                                    "cutoff_style": cutoff_style,
                                                    "r_cut": r_cuts[
                                                        molecule
                                                    ].to_value("nm"),
                                                }
                                                total_statepoints.append(
                                                    statepoint
                                                )

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
