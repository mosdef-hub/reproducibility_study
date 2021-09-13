"""Script for analysing the results from the 16 different seeds from the 11 systems. It also parses the hdf5 format trajectory stored in each output analysis folder (obtained by executing conv_traj.py before this script) to get the RDFs."""
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
    data_path = "analysis_data"
    os.chdir(data_path)

    project = signac.get_project()

    for (molecule, ensemble, temperature, pressure), group in project.groupby(
        (
            "molecule",
            "ensemble",
            "temperature",
            "pressure",
        )
    ):
        print(molecule, ensemble, temperature, pressure)
        os.chdir(
            "{}_{}_{}K_{}kPa".format(molecule, ensemble, temperature, pressure)
        )
        base_dir = os.getcwd()
        if ensemble == "NPT" and molecule == "methaneUA":

            density_list = []
            for job in group:
                os.chdir(job.ws)
                prod_run_files = sorted(glob("run*prod*"))
                if len(prod_run_files) < 4:
                    print(
                        "warning, only {}  prod cycles complete for {} {} {} {} {}".format(
                            len(prod_run_files),
                            job,
                            molecule,
                            ensemble,
                            temperature,
                            pressure,
                        )
                    )
                density_list.append(avg_one_seed_density(prod_run_files))
            print(len(traj_list))
            filtered_density_list = list(filter(None, density_list))
            output_string = "The average density is {} g/ml with SEM {} from {} samples".format(
                np.mean(filtered_density_list),
                np.std(filtered_density_list)
                / np.sqrt(len(filtered_density_list)),
                len(filtered_density_list),
            )
            print(output_string)
            os.chdir(base_dir)
            text_file = open("density_results.txt", "w")
            text_file.write(output_string)
            text_file.close()

            comb_traj = md.load_hdf5("comb_traj.h5")
            print(comb_traj)
        elif ensemble == "GEMC-NVT" and molecule == "methaneUA":
            for job in group:
                print(job)

        os.chdir("..")


def avg_one_seed_density(prod_run_files):
    """For a one particular seed, read all the prod run files and provide the average density (g/ml) for one seed."""
    if len(prod_run_files) == 0:
        return None
    os.system(
        "grep 'specific density                        ' run.prod* | awk '{print $5'} > temp_density.txt"
    )
    density_list_one_seed = np.genfromtxt("temp_density.txt")
    # now delete the temp file
    os.remove("temp_density.txt")
    return np.mean(density_list_one_seed)


if __name__ == "__main__":
    main()
