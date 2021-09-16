import signac
import os
import ele
from unyt import g, mole, angstrom
import numpy as np

root = os.path.join(os.getcwd(), "../../../src")
project = signac.get_project(root)

results = {}

NA = 6.032E23

C = ele.element_from_symbol("C")
H = ele.element_from_symbol("H")

mw = C.mass + H.mass * 4
mw *= g / mole

density = np.array([])

for job in project.find_jobs({"molecule": "methaneUA", "engine": "cassandra", "ensemble": "NPT"}):
    results[job.statepoint()["replica"]] = {}
    density = np.append(density, job.document["mean_density_box1"])
