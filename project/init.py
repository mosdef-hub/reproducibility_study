"""Initialize signac statepoints."""
import itertools
import os

import signac


def dict_product(dd):
    """Return the product of the key/values of a dictionary."""
    keys = dd.keys()
    for element in itertools.product(*dd.values()):
        yield dict(zip(keys, element))


simulation_engines = [
    "cassandra",
    "mcccs",
    "gomc",
    "gromacs",
    "hoomd",
    "lammps",
]
forcefields = ["Trappe_UA", "spce", "oplsaa"]
densities = [1.0]  # units
temperatures = [1.0]  # units
molecules = ["methaneUA", "pentaneUA", "benzeneUA", "waterSPC/E", "ethanolAA"]
N_molecules = [1000]
ensembles = ["NPT", "GEMC-NVT"]
replicas = range(5)
params = {
    "simulation_engine": simulation_engines,
    "forcefield": forcefields,
    "density": densities,
    "temperature": temperatures,
    "molecule": molecules,
    "N_molecules": N_molecules,
    "production_ensemble": ensembles,
    "replica": replicas,
}

pr_root = os.path.join(os.getcwd(), "src")
pr = signac.get_project(pr_root)

# ignore statepoints that are not being tested (gemc only for methane, pentane)
# filter the list of dictionaries
total_statepoints = list()
indices_to_remove = set()
for sp in dict_product(params):
    total_statepoints.append(sp)

for i, sp in enumerate(total_statepoints):
    # filter gemc ensembles from md engines
    if sp["production_ensemble"] == "GEMC-NVT" and sp["simulation_engine"] in [
        "gromacs",
        "hoomd",
        "lammps",
    ]:
        indices_to_remove.add(i)

    # filter gemc production_ensembles from all but 2 molecules
    if sp["production_ensemble"] == "GEMC-NVT" and sp["molecule"] not in [
        "methaneUA",
        "pentaneUA",
    ]:
        indices_to_remove.add(i)

    # filter UA ffs from AA/water models
    if sp["forcefield"] == "Trappe_UA" and sp["molecule"] in [
        "waterSPC/E",
        "ethanolAA",
    ]:
        indices_to_remove.add(i)

    # filter oplsaa FFs from UA and spc/e water
    if sp["forcefield"] == "oplsaa" and sp["molecule"] not in ["ethanolAA"]:
        indices_to_remove.add(i)

    # filter non-water molecules from spce ff
    if sp["forcefield"] == "spce" and sp["molecule"] != "waterSPC/E":
        indices_to_remove.add(i)

# now reverse sort the set and remove from inital list
# must be reverse sorted to remove indices on the list in place
# otherwise the list will change size and the indices would change
sorted_indicies_to_delete = sorted(list(indices_to_remove), reverse=True)
for idx in sorted_indicies_to_delete:
    del total_statepoints[idx]

for sp in total_statepoints:
    pr.open_job(
        statepoint=sp,
    ).init()
