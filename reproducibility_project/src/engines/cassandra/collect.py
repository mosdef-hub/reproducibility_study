"""Gather information for post simulation analysis across statepoints."""

import os

import numpy as np
import signac

root = os.path.join(os.getcwd(), "../../../src")
project = signac.get_project(root)

results = {}

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
