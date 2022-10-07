"""Script for finding the job ids of specific simulations."""
# It also parses the gsd format trajectory stored in each output analysis folder (obtained by executing conv_traj.py before this script) to get the RDFs."""
import os
import shutil
from glob import glob

import freud
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import signac
from scipy import stats


def main():

    project = signac.get_project()

    for (
        molecule,
        ensemble,
        temperature,
        pressure,
        cutoff_style,
        long_range_correction,
    ), group in project.groupby(
        (
            "molecule",
            "ensemble",
            "temperature",
            "pressure",
            "cutoff_style",
            "long_range_correction",
        )
    ):

        print("-----------------------------------------------------")

        if molecule == "methaneUA" and ensemble == "NPT":
            print(
            molecule,
            ensemble,
            temperature,
            pressure,
            cutoff_style,
            long_range_correction,
            )
            
            for job in group:
                print(job.sp.replica, job)



if __name__ == "__main__":
    main()
