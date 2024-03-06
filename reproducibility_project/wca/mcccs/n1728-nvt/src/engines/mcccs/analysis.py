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
    data_path = "analysis_data"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    # delete the folder manually

    os.chdir(data_path)

    project = signac.get_project()

    for (
        engine,
        molecule,
        ensemble,
        temperature,
        pressure,
        cutoff_style,
        long_range_correction,
    ), group in project.groupby(
        (
            "engine",
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
            engine,
            molecule,
            ensemble,
            temperature,
            pressure,
            cutoff_style,
            long_range_correction,
        )
        if not os.path.isdir(
            "{}_{}_{}_{}K_{}kPa_cutoff_{}_lrc_{}".format(
                engine,
                molecule,
                ensemble,
                temperature,
                pressure,
                cutoff_style,
                str(long_range_correction),
            )
        ):
            os.makedirs(
                "{}_{}_{}_{}K_{}kPa_cutoff_{}_lrc_{}".format(
                    engine,
                    molecule,
                    ensemble,
                    temperature,
                    pressure,
                    cutoff_style,
                    str(long_range_correction),
                )
            )
        os.chdir(
            "{}_{}_{}_{}K_{}kPa_cutoff_{}_lrc_{}".format(
                engine,
                molecule,
                ensemble,
                temperature,
                pressure,
                cutoff_style,
                str(long_range_correction),
            )
        )

        base_dir = os.getcwd()
        if "NVT" in ensemble and engine == "mcccs":
            density_list = []
            for job in group:
                os.chdir(job.ws)
                # print(job)
                prod_run_files = sorted(glob("run*prod*"))
                if len(prod_run_files) < 4:
                    print(
                        "warning, {} engine only {}  prod cycles complete for {} {} {} {} {} {} {}".format(
                            engine,
                            len(prod_run_files),
                            job,
                            molecule,
                            ensemble,
                            temperature,
                            pressure,
                            cutoff_style,
                            long_range_correction,
                        )
                    )
                density_list.append(avg_one_seed_density_box1(prod_run_files))

                # print(density_list)
            filtered_density_list = list(filter(None, density_list))
            output_string = "The average density is {} number/Ang^3 with SEM {} from {} samples".format(
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

        elif ensemble == "GEMC-NVT":  # and molecule == "methaneUA":
            density_list_box1 = []
            density_list_box2 = []
            for job in group:
                # print("This is GEMC ",job)
                os.chdir(job.ws)
                prod_run_files = sorted(glob("run*prod*"))
                if len(prod_run_files) < 4:
                    print(
                        "warning, only {}  prod cycles complete for {} {} {} {} {} {} {}".format(
                            len(prod_run_files),
                            job,
                            molecule,
                            ensemble,
                            temperature,
                            pressure,
                            cutoff_style,
                            long_range_correction,
                        )
                    )
                density_list_box1.append(
                    avg_one_seed_density_box1(prod_run_files)
                )
                density_list_box2.append(
                    avg_one_seed_density_box2(prod_run_files)
                )

                # print(density_list)
            filtered_density_list_box1 = list(filter(None, density_list_box1))
            filtered_density_list_box2 = list(filter(None, density_list_box2))

            output_string_box1 = "Box1 The average density is {} number/Ang^3 with SEM {} from {} samples".format(
                np.mean(filtered_density_list_box1),
                np.std(filtered_density_list_box1)
                / np.sqrt(len(filtered_density_list_box1)),
                len(filtered_density_list_box1),
            )
            output_string_box2 = "Box2 The average density is {} number/Ang^3 with SEM {} from {} samples".format(
                np.mean(filtered_density_list_box2),
                np.std(filtered_density_list_box2)
                / np.sqrt(len(filtered_density_list_box2)),
                len(filtered_density_list_box2),
            )
            print(output_string_box1)
            print(output_string_box2)
            os.chdir(base_dir)
            with open("density_results.txt", "w") as text_file:
                text_file.write(output_string_box1 + "\n" + output_string_box2)
            text_file.close()

        os.chdir("..")


def avg_one_seed_density_box1(prod_run_files):
    """For a one particular seed, read all the prod run files and provide the average density (number/Ang^3) for one seed."""
    if len(prod_run_files) == 0:
        return None
    os.system(
        "grep 'number density of type    1 PEC1' run.prod* | awk '{print $10}' > temp_density.txt"
    )
    # os.system("grep 'specific density                        ' run.prod* | awk '{print $6}' ")
    if os.path.exists("temp_density.txt"):
        try:
            density_list_one_seed = np.genfromtxt("temp_density.txt") / 1000
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp_density.txt")
        return np.mean(density_list_one_seed)
    else:
        return None


def avg_one_seed_density_box2(prod_run_files):
    """For a one particular seed, read all the prod run files and provide the average density (number/Ang^3) for one seed."""
    if len(prod_run_files) == 0:
        return None
    os.system(
        "grep 'number density of type    1 PEC1' run.prod* | awk '{print $10}' > temp_density.txt"
    )
    # os.system("grep 'specific density                        ' run.prod* | awk '{print $6}' ")
    if os.path.exists("temp_density.txt"):
        try:
            density_list_one_seed = np.genfromtxt("temp_density.txt") / 1000
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp_density.txt")
        return np.mean(density_list_one_seed)
    else:
        return None


if __name__ == "__main__":
    main()
