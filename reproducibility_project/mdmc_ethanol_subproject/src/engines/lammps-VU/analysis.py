"""Script for analysing the results from the 16 different seeds from the systems."""
# It also parses the gsd format trajectory stored in each output analysis folder (obtained by executing conv_traj.py before this script) to get the RDFs."""
import os
import shutil
from glob import glob

import freud
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
import signac
from scipy import stats

from reproducibility_project.src.analysis.equilibration import is_equilibrated


def main():
    """Read run files from all the independent seeds and get the final density and RDFs."""
    data_path = "analysis_data"
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
        print(
            molecule,
            ensemble,
            temperature,
            pressure,
            cutoff_style,
            long_range_correction,
        )
        if not os.path.isdir(
            "{}_{}_{}K_{}kPa_cutoff_{}_lrc_{}".format(
                molecule,
                ensemble,
                temperature,
                pressure,
                cutoff_style,
                str(long_range_correction),
            )
        ):

            os.makedirs(
                "{}_{}_{}K_{}kPa_cutoff_{}_lrc_{}".format(
                    molecule,
                    ensemble,
                    temperature,
                    pressure,
                    cutoff_style,
                    str(long_range_correction),
                )
            )
        os.chdir(
            "{}_{}_{}K_{}kPa_cutoff_{}_lrc_{}".format(
                molecule,
                ensemble,
                temperature,
                pressure,
                cutoff_style,
                str(long_range_correction),
            )
        )

        base_dir = os.getcwd()
        if "NPT" in ensemble:
            density_list = []
            for job in group:
                os.chdir(job.ws)
                # print(job)
                density_list.append(avg_one_seed_density_box1(job))

            filtered_density_list = list(
                filter(lambda x: x is not None, density_list)
            )
            output_string = "The average density is {} g/ml with SEM {} from {} samples".format(
                np.mean(filtered_density_list),
                np.std(filtered_density_list)
                / np.sqrt(len(filtered_density_list)),
                len(filtered_density_list),
            )
            print(output_string)
            os.chdir(base_dir)
            with open("density_results.txt", "w") as text_file:
                text_file.write(output_string)
            text_file.close()

        os.chdir("..")


def avg_one_seed_density_box1(job):
    """For a one particular seed, read all the prod run files and provide the average density (g/ml) for one seed."""
    try:
        data = pd.read_csv(job.ws + "/log-npt.txt", delimiter=" ", header=0)
    except:
        return None
    attr = "density"
    passed, start, step, neff = is_equilibrated(data[attr], 0.2, 1)
    if passed:
        temp = job.sp.temperature
        thresh_frac = round((len(data) - start) / len(data), 3)
        started_data = data[attr][start:]
        # uncorr_data = started_data[0::math.ceil(step)]
        uncorr_data = started_data[0::]

        return np.mean(uncorr_data)


if __name__ == "__main__":
    main()
