"""Script for deleting certain files in a state point folder."""
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
    """Load the project and finding the simulation folder."""
    project = signac.get_project()

    for (
        molecule,
        ensemble,
        temperature,
        pressure,
        cutoff_style,
        long_range_correction,
        engine,
    ), group in project.groupby(
        (
            "molecule",
            "ensemble",
            "temperature",
            "pressure",
            "cutoff_style",
            "long_range_correction",
            "engine",
        )
    ):
        print("-----------------------------------------------------")

        if (
            molecule == "methaneUA"
            and ensemble == "NPT-gigantic"
            and engine == "mcccs"
        ):
            print(
                molecule,
                ensemble,
                temperature,
                pressure,
                cutoff_style,
                long_range_correction,
                engine,
            )

            for job in group:
                os.chdir(job.workspace())
                # print(os.listdir())
                # os.system("cat sig*state*json")
                os.system("rm -rf *prod*")
                if job.doc.prod_replicates_done > 0:
                    job.doc.prod_replicates_done = 0
                # print(os.listdir())


if __name__ == "__main__":
    main()
