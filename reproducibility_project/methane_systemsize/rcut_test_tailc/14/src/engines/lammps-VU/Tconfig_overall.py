"""Script for analysing the results from the 16 different seeds from the 11 systems."""
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
    """Read run files from all the independent seeds and get the final density and RDFs."""
    data_path = "Tconfig_analysis_data"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    # delete the folder manually

    os.chdir(data_path)

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
        print(
            molecule,
            ensemble,
            temperature,
            pressure,
            cutoff_style,
            long_range_correction,
            engine,
        )
        if not os.path.isdir(
            "{}_{}_{}K_{}kPa_cutoff_{}_lrc_{}_{}".format(
                molecule,
                ensemble,
                temperature,
                pressure,
                cutoff_style,
                str(long_range_correction),
                engine,
            )
        ):
            os.makedirs(
                "{}_{}_{}K_{}kPa_cutoff_{}_lrc_{}_{}".format(
                    molecule,
                    ensemble,
                    temperature,
                    pressure,
                    cutoff_style,
                    str(long_range_correction),
                    engine,
                )
            )
        os.chdir(
            "{}_{}_{}K_{}kPa_cutoff_{}_lrc_{}_{}".format(
                molecule,
                ensemble,
                temperature,
                pressure,
                cutoff_style,
                str(long_range_correction),
                engine,
            )
        )

        base_dir = os.getcwd()

        if engine == "lammps-VU":
            Tconfig_list = []
            for job in group:
                os.chdir(job.ws)
                Tconfig = np.genfromtxt("Tconfig.txt")
                Tconfig_list.append(Tconfig)

            os.chdir(base_dir)
            avg_Tconfig = np.mean(Tconfig_list, axis=0)
            sem_Tconfig = (
                np.std(Tconfig_list, axis=0) / (len(Tconfig_list)) ** 0.5
            )

            np.savetxt("Tconfig.txt", np.vstack((avg_Tconfig, sem_Tconfig)))

        os.chdir("..")


if __name__ == "__main__":
    main()
