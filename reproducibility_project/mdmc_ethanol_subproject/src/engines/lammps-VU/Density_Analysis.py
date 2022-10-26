"""Quick Analysis to output average densities locally for LAMMPS-VU."""
import math
import sys

import matplotlib.pyplot as plt
import mbuild as mb
import mdtraj as md
import numpy as np
import pandas as pd
import pymbar
import signac

# sys.path.append('reproducibility_study-1/reproducibility_project/src/analysis/')
from reproducibility_project.src.analysis.equilibration import is_equilibrated

if __name__ == "__main__":
    path = "./"
    project = signac.get_project(path)
    print(project.id)
    analysis_jobs_list = []
    for job in project:
        if job.sp.engine == "lammps-VU":
            if job.sp.molecule == "ethanolAA":
                if job.isfile("log-npt.txt"):
                    analysis_jobs_list.append(job)

    print(len(analysis_jobs_list))
    min_thresh = 1
    means = {280.0: [], 300.0: [], 320.0: [], 400.0: []}
    stds = {280.0: [], 300.0: [], 320.0: [], 400.0: []}
    for job in analysis_jobs_list:
        data = pd.read_csv(job.ws + "/log-npt.txt", delimiter=" ", header=0)
        attr = "density"
        passed, start, step, neff = is_equilibrated(data[attr], 0.2, 1)
        if passed:
            temp = job.sp.temperature
            print(f"{job.sp.molecule} {job.id}\n")
            print(
                "Step size is: {}\nStarting Frame: {}\nNumber of Effective Frames: {}".format(
                    step, start, neff
                )
            )
            thresh_frac = round((len(data) - start) / len(data), 3)
            print("Threshold Fraction is {}".format(thresh_frac))
            print(np.mean(data[attr]), np.std(data[attr]))
            started_data = data[attr][start:]
            print(
                "Length of usable data ",
                len(started_data[0 :: math.ceil(step)]),
            )
            print("\n______________________________\n")
            if thresh_frac < min_thresh:
                min_thresh = thresh_frac
                thresh_id = job.id
            # uncorr_data = started_data[0::math.ceil(step)]
            uncorr_data = started_data[0::]
            means[temp].append(np.mean(uncorr_data))
            stds[temp].append(np.std(uncorr_data))

    for temp in means:
        if len(means[temp]) >= 2:
            print(
                "Density is: ",
                np.mean(means[temp]),
                np.std(means[temp]),
                f" at T={temp}",
            )
    print("Minimum Threshhold is: {} id: {}".format(min_thresh, thresh_id))
    # fig = plt.hist(means, bins=9)
    # plt.savefig('Benzene_Histogram.pdf')
    # plt.scatter(np.arange(0,len(means)), np.array(means))
    # plt.savefig('Means of Pentane')
