"""Initialize signac statepoints."""
import itertools
import os

import signac


def dict_product(dd):
    """Create a dictionary that is the product of the key,value of a dict."""
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
forcefields = ["Trappe_UA"]
densities = [1.0]  # units
temperatures = [1.0]  # units
molecules = ["methane", "ethane", "propane", "butane", "pentane", "hexane"]
N_molecules = [1000]
replicas = range(5)
params = {
    "simulation_engine": simulation_engines,
    "forcefield": forcefields,
    "density": densities,
    "temperature": temperatures,
    "molecule": molecules,
    "N_molecules": N_molecules,
    "replica": replicas,
}

pr_root = os.path.join(os.getcwd(), "src")
pr = signac.get_project(pr_root)
for sp in dict_product(params):
    job = pr.open_job(sp).init()
