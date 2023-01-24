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


def Power_Analysis(means_S1, means_S2, delta, zbeta=1.645, zalpha=1.960):
    """Generate number of samples for required power of study."""
    return (
        (np.std(means_S1) ** 2 + np.std(means_S2) ** 2)
        * (zalpha + zbeta) ** 2
        / delta**2
    )


if __name__ == "__main__":
    path = "./"
    project = signac.get_project(path)
    print(project.id)
    analysis_jobs_list = []
    for job in project:
        if job.sp.engine == "lammps-VU":
            if job.sp.molecule == "methaneUA":
                if job.isfile("log-npt.txt"):
                    analysis_jobs_list.append(job)

    print(len(analysis_jobs_list))
    min_thresh = 1
    N_liq_molecules = [450, 600, 900, 1800]
    means = {}
    stds = {}
    [
        (means.__setitem__(key, []), stds.__setitem__(key, []))
        for key in N_liq_molecules
    ]
    for job in analysis_jobs_list:
        data = pd.read_csv(job.ws + "/log-npt.txt", delimiter=" ", header=0)
        attr = "density"
        passed, start, step, neff = is_equilibrated(data[attr], 0.2, 1)
        if passed:
            mols = job.sp.N_liquid
            print(f"{mols} {job.id}\n")
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
            means[mols].append(np.mean(uncorr_data))
            stds[mols].append(np.std(uncorr_data))
        else:
            print(f"job {job.id} is not equilibrated")

    for mols in means:
        if len(means[mols]) >= 2:
            print(
                "Density is: ",
                np.mean(means[mols]),
                np.std(means[mols]),
                f" at Nmols={mols}",
            )
    # fig = plt.hist(means, bins=9)
    # plt.savefig('Benzene_Histogram.pdf')
    # plt.scatter(np.arange(0,len(means)), np.array(means))
    # plt.savefig('Means of Pentane')

    # Calcuate number of samples necessary to see differences between means
    delta = 4e-5
    print(
        f"Number of samples needed to see a difference of {delta} is "
        f"{Power_Analysis(means[600], means[900], delta=delta, zalpha=1.645, zbeta=1.282)}"
    )
