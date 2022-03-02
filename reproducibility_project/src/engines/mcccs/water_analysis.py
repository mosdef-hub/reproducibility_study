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
    data_path = "water_data"
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
        if molecule != "waterSPCE":
            continue

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
        if ensemble == "NPT":
            density_list = []
            cp_list = []
            press_list = []  # list of pressure from run.prod files
            ener_list = []  # list of energy from run.prod fles
            len_list = []  # list of length from run.prod files
            energies_list = []  # reading vectors of energies for 16 seeds
            length_list = []  # reading vectors of lengths for 16 seeds
            lengths = []  # averages from fort12 files
            pressures = []  # averages from fort12 files
            energies = []  # averages from fort12 files
            trans_max_disp = []
            rot_max_disp = []
            vol_max_disp = []
            for job in group:
                os.chdir(job.ws)
                # print(job)
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
                density_list.append(avg_one_seed_density_box1(prod_run_files))
                cp_list.append(avg_one_seed_cp_box1(prod_run_files))
                press_list.append(avg_one_seed_press_box1(prod_run_files))
                ener_list.append(avg_one_seed_ener_box1(prod_run_files))
                len_list.append(avg_one_seed_len_box1(prod_run_files))
                fort12data = np.genfromtxt("log-npt.txt", names=True)
                # print("The shape is {}".format(
                energies_list.append(
                    np.reshape(fort12data["potential_energy"], (120000, 1))
                )
                length_list.append(np.reshape(fort12data["a"], (120000, 1)))
                lengths.append(
                    np.mean(np.reshape(fort12data["a"], (120000, 1)))
                )
                pressures.append(
                    np.mean(np.reshape(fort12data["pressure"], (120000, 1)))
                )
                energies.append(
                    np.mean(
                        np.reshape(fort12data["potential_energy"], (120000, 1))
                    )
                )
                trans_max_disp.append(get_trans_max_disp())
                rot_max_disp.append(get_rot_max_disp())
                vol_max_disp.append(get_vol_max_disp())
                # energies_list.append(fort12data['potential_energy'])

                # print(density_list)
            print(trans_max_disp)
            filtered_density_list = list(filter(None, density_list))
            filtered_cp_list = list(filter(None, cp_list))
            filtered_press_list = list(filter(None, press_list))
            filtered_ener_list = list(filter(None, ener_list))
            filtered_len_list = list(filter(None, len_list))
            output_string = "The average density from run files is {} g/ml with SEM {} from {} samples".format(
                np.mean(filtered_density_list),
                np.std(filtered_density_list)
                / np.sqrt(len(filtered_density_list)),
                len(filtered_density_list),
            )
            output_string += "\n"
            output_string += "The average cp from run files is {} J/Kmol with SEM {} from {} samples".format(
                np.mean(filtered_cp_list),
                np.std(filtered_cp_list) / np.sqrt(len(filtered_cp_list)),
                len(filtered_cp_list),
            )
            output_string += "\n"

            output_string += "The average energy from run files is {} K with SEM {} from {} samples".format(
                np.mean(filtered_ener_list),
                np.std(filtered_ener_list) / np.sqrt(len(filtered_ener_list)),
                len(filtered_ener_list),
            )
            output_string += "\n"

            output_string += "The average energy from run files is {} kJ/mol with SEM {} from {} samples".format(
                0.008314410016255 * np.mean(filtered_ener_list),
                0.008314410016255
                * np.std(filtered_ener_list)
                / np.sqrt(len(filtered_ener_list)),
                len(filtered_ener_list),
            )
            output_string += "\n"

            output_string += "The average press from run files is {} kPa with SEM {} from {} samples".format(
                np.mean(filtered_press_list),
                np.std(filtered_press_list) / np.sqrt(len(filtered_press_list)),
                len(filtered_press_list),
            )
            output_string += "\n"

            output_string += "The average len from run files is {} nm with SEM {} from {} samples".format(
                np.mean(filtered_len_list) / 10,
                np.std(filtered_len_list)
                / 10
                / np.sqrt(len(filtered_len_list)),
                len(filtered_len_list),
            )
            output_string += "\n"

            output_string += "The average box length from fort12 files is {} nm with SEM {} from {} samples".format(
                np.mean(lengths),
                np.std(lengths) / np.sqrt(len(lengths)),
                len(lengths),
            )

            output_string += "\n"
            output_string += "The average box pressure from fort12 files is {} kPa with SEM {} from {} samples".format(
                np.mean(pressures),
                np.std(pressures) / np.sqrt(len(pressures)),
                len(pressures),
            )

            output_string += "\n"
            output_string += "The average box energy from fort12 files is {} kJ/mol with SEM {} from {} samples".format(
                np.mean(energies),
                np.std(energies) / np.sqrt(len(energies)),
                len(energies),
            )
            output_string += "\n"
            output_string += "The average trans max disp is {}  with SEM {} from {} samples".format(
                np.mean(trans_max_disp),
                np.std(trans_max_disp) / np.sqrt(len(trans_max_disp)),
                len(trans_max_disp),
            )

            output_string += "\n"
            output_string += "The average rot  max disp is {}  with SEM {} from {} samples".format(
                np.mean(rot_max_disp),
                np.std(rot_max_disp) / np.sqrt(len(rot_max_disp)),
                len(rot_max_disp),
            )

            output_string += "\n"
            output_string += "The average vol max disp is {}  with SEM {} from {} samples".format(
                np.mean(vol_max_disp),
                np.std(vol_max_disp) / np.sqrt(len(vol_max_disp)),
                len(vol_max_disp),
            )

            print(output_string)
            os.chdir(base_dir)
            with open("density_results.txt", "w") as text_file:
                text_file.write(output_string)

            text_file.close()
            textfile = open("trans_max_disp.txt", "w")
            for element in trans_max_disp:
                textfile.write(str(element) + "\n")
            textfile.close()
            plt.plot(trans_max_disp)
            plt.xlabel("Seed")
            plt.ylabel("Max displacement (\AA)")
            plt.tight_layout()
            plt.savefig("Trans max displacement_{}.png".format(temperature))
            plt.close()

            textfile = open("rot_max_disp.txt", "w")
            for element in rot_max_disp:
                textfile.write(str(element) + "\n")
            textfile.close()
            plt.plot(rot_max_disp)
            plt.xlabel("Seed")
            plt.ylabel("Max displacement (radians)")
            plt.tight_layout()
            plt.savefig("Rot max displacement_{}.png".format(temperature))
            plt.close()

            textfile = open("vol_max_disp.txt", "w")
            for element in vol_max_disp:
                textfile.write(str(element) + "\n")
            textfile.close()
            plt.plot(vol_max_disp)
            plt.xlabel("Seed")
            plt.ylabel("Max displacement (\AA^3)")
            plt.tight_layout()
            plt.savefig("Vol max displacement_{}.png".format(temperature))

            plt.close()
            print(len(energies_list))
            energies_list = np.concatenate(energies_list, axis=1)
            length_list = np.concatenate(length_list, axis=1)
            plt.plot(energies_list)
            plt.xlabel("MC cycles")
            plt.ylabel("Energy (kJ/mol)")
            plt.tight_layout()
            plt.savefig("water_energy_evolution_{}.png".format(temperature))
            plt.close()
            plt.plot(length_list)
            plt.xlabel("MC cycles")
            plt.ylabel("Box length (nm)")
            plt.tight_layout()
            plt.savefig("water_box_length_evolution_{}.png".format(temperature))
            plt.close()

        os.chdir("..")


def avg_one_seed_density_box1(prod_run_files):
    """For a one particular seed, read all the prod run files and provide the average density (g/ml) for one seed."""
    if len(prod_run_files) == 0:
        return None
    os.system(
        "grep 'specific density                        ' run.prod* | awk '{print $6}' > temp_density.txt"
    )
    # os.system("grep 'specific density                        ' run.prod* | awk '{print $6}' ")
    if os.path.exists("temp_density.txt"):
        try:
            density_list_one_seed = np.genfromtxt("temp_density.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp_density.txt")
        return np.mean(density_list_one_seed)
    else:
        return None


def avg_one_seed_cp_box1(prod_run_files):
    """For a one particular seed, read all the prod run files and provide the average cp (j/Kmol) for one seed."""
    if len(prod_run_files) == 0:
        return None
    os.system(
        "grep 'Cp residual(J/Kmol)' run.prod* | awk '{print $5}' > temp_cp.txt"
    )
    # os.system("grep 'specific density                        ' run.prod* | awk '{print $6}' ")
    if os.path.exists("temp_cp.txt"):
        try:
            cp_list_one_seed = np.genfromtxt("temp_cp.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp_cp.txt")
        return np.mean(cp_list_one_seed)
    else:
        return None


def avg_one_seed_press_box1(prod_run_files):
    """For a one particular seed, read all the prod run files and provide the average pressure (kPa) for one seed."""
    if len(prod_run_files) == 0:
        return None
    os.system(
        "grep 'pressure                                      ' run.prod* | awk '{print $5}' > temp_press.txt"
    )
    # os.system("grep 'specific density                        ' run.prod* | awk '{print $6}' ")
    if os.path.exists("temp_press.txt"):
        try:
            press_list_one_seed = np.genfromtxt("temp_press.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp_press.txt")
        return np.mean(press_list_one_seed)
    else:
        return None


def avg_one_seed_ener_box1(prod_run_files):
    """For a one particular seed, read all the prod run files and provide the average energy (K) for one seed."""
    if len(prod_run_files) == 0:
        return None
    os.system(
        "grep 'Total energy \[K per system ' run.prod* | awk '{print $12}' > temp_ener.txt"
    )
    # os.system("grep 'specific density                        ' run.prod* | awk '{print $6}' ")
    if os.path.exists("temp_ener.txt"):
        try:
            ener_list_one_seed = np.genfromtxt("temp_ener.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp_ener.txt")
        return np.mean(ener_list_one_seed)
    else:
        return None


def avg_one_seed_len_box1(prod_run_files):
    """For a one particular seed, read all the prod run files and provide the average len  (nm) for one seed."""
    if len(prod_run_files) == 0:
        return None
    os.system(
        "grep 'boxlength                                       \[' run.prod* | awk '{print $5}' > temp_len.txt"
    )
    # os.system("grep 'specific density                        ' run.prod* | awk '{print $6}' ")
    if os.path.exists("temp_len.txt"):
        try:
            len_list_one_seed = np.genfromtxt("temp_len.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp_len.txt")
        return np.mean(len_list_one_seed)
    else:
        return None


def avg_one_seed_density_box2(prod_run_files):
    """For a one particular seed, read all the prod run files and provide the average density (g/ml) for one seed."""
    if len(prod_run_files) == 0:
        return None
    os.system(
        "grep 'specific density                        ' run.prod* | awk '{print $7}' > temp_density.txt"
    )
    # os.system("grep 'specific density                        ' run.prod* | awk '{print $6}' ")
    if os.path.exists("temp_density.txt"):
        try:
            density_list_one_seed = np.genfromtxt("temp_density.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp_density.txt")
        return np.mean(density_list_one_seed)
    else:
        return None


def get_trans_max_disp():
    """Get the max disp for translation moves."""
    f = open("config1a.dat.prod0", "r")
    k = 0
    for line in f:
        if k == 2:
            return float(line.split()[0])
        k += 1


def get_rot_max_disp():
    """Get the max disp for rotation moves."""
    f = open("config1a.dat.prod0", "r")
    k = 0
    for line in f:
        if k == 3:
            return float(line.split()[0])
        k += 1


def get_vol_max_disp():
    """Get the max disp for vol moves."""
    f = open("config1a.dat.prod0", "r")
    k = 0
    for line in f:
        if k == 5:
            return float(line.split()[0])
        k += 1


if __name__ == "__main__":
    main()
