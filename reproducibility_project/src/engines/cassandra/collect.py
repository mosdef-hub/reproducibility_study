"""Gather information for post simulation analysis across statepoints."""

import os

import ele
import numpy as np
import signac
from unyt import angstrom, g, mole

root = os.path.join(os.getcwd(), "../../../src")
project = signac.get_project(root)

results = {}

NA = 6.032e23

C = ele.element_from_symbol("C")
H = ele.element_from_symbol("H")

mw = C.mass + H.mass * 4
mw *= g / mole

density = np.array([])

for job in project.find_jobs(
    {"molecule": "methaneUA", "engine": "cassandra", "ensemble": "NPT"}
):
    results[job.statepoint()["replica"]] = {}
    density = np.append(density, job.document["mean_density_box1"])


mean = density.mean()
stdev = density.std()

print(density)
print(mean)
print(stdev)
