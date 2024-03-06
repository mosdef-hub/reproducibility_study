"""Script for analysing the results from the 16 different seeds from the 11 systems."""

# It also parses the gsd format trajectory stored in each output analysis folder (obtained by executing conv_traj.py before this script) to get the RDFs."""
import os
import shutil
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import signac
from scipy import stats


def main():
    """Read run files from all the independent seeds and get the final density and RDFs."""
    data_path = "spe_data"
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
        if ensemble == "NPT":
            density_list = []
            for job in group:
                if job.sp.engine != "mcccs":
                    continue
                os.chdir(job.ws)
                potential_energy = get_pe()
                lj_inter_energy = get_lj_inter()
                lj_intra_energy = get_lj_intra()
                tail_correction = get_tailc()
                coulomb_energy = get_coul()
                bond_energy = get_bond()
                angle_energy = get_angle()
                dihedral_energy = get_torsion()
                print("working on ", job)
                csv_pe = potential_energy
                csv_lj_energy = lj_inter_energy + lj_intra_energy
                csv_tail_energy = tail_correction
                csv_coulomb_energy = float("NaN")
                csv_kspace_energy = float("NaN")
                # csv_kspace_energy = 0
                csv_total_coulombic = coulomb_energy
                csv_pair_energy = csv_lj_energy + csv_total_coulombic
                csv_bond_energy = bond_energy
                csv_angle_energy = angle_energy
                csv_dihedral_energy = dihedral_energy
                csv_mol_energy = bond_energy + angle_energy + dihedral_energy
                csv_intramolecular_energy = float("NaN")
                csv_intermolecular_energy = float("NaN")

                # csv_intermolecular_energy = lj_inter_energy
                # csv_intramolecular_energy = lj_intra_energy
                # csv_short_range = lj_inter_energy + lj_intra_energy

                output_string = ""

                print(output_string)
                energies = 0.008314462618 * np.array(
                    [
                        [
                            csv_pe,
                            csv_lj_energy,
                            csv_tail_energy,
                            csv_total_coulombic,
                            csv_coulomb_energy,
                            csv_kspace_energy,
                            csv_pair_energy,
                            csv_bond_energy,
                            csv_angle_energy,
                            csv_dihedral_energy,
                            csv_mol_energy,
                            csv_intramolecular_energy,
                            csv_intermolecular_energy,
                        ]
                    ]
                )
                header = "potential_energy \t tot_vdw_energy \t tail_energy \t tot_electrostatics \t short_range_electrostatics \t long_range_electrostatics \t tot_pair_energy \t bonds_energy \t angles_energy \t dihedrals_energy \t tot_bonded_energy \t intramolecular_energy \t intermolecular_energy"

                df = pd.DataFrame(
                    energies,
                    columns=[
                        "potential_energy",
                        "tot_vdw_energy",
                        "tail_energy",
                        "tot_electrostatics",
                        "short_range_electrostatics",
                        "long_range_electrostatics",
                        "tot_pair_energy",
                        "bonds_energy",
                        "angles_energy",
                        "dihedrals_energy",
                        "tot_bonded_energy",
                        "intramolecular_energy",
                        "intermolecular_energy",
                    ],
                )
                df.to_csv("log-spe.txt", header=True, index=False, sep=",")

                os.chdir(base_dir)
                np.savetxt(
                    "log-spe.txt", energies, header=header, delimiter="\t"
                )

        os.chdir("..")


def get_pe():
    """Get PE for the system."""
    os.system("grep 'total energy' run.melt | awk '{print $3}' > temp_pe.txt")
    if os.path.exists("temp_pe.txt"):
        try:
            pe_one_seed = np.genfromtxt("temp_pe.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp_pe.txt")
        return np.mean(pe_one_seed[0])
    else:
        return None


def get_lj_inter():
    """Get inter LJ energy for the system."""
    os.system(
        "grep 'inter lj energy' run.melt | awk '{print $4}' > temp_lj_inter.txt"
    )
    if os.path.exists("temp_lj_inter.txt"):
        try:
            lj_inter_one_seed = np.genfromtxt("temp_lj_inter.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp_lj_inter.txt")
        return np.mean(lj_inter_one_seed[0])
    else:
        return None


def get_lj_intra():
    """Get intra LJ energy for the system."""
    os.system(
        "grep 'intra lj energy' run.melt | awk '{print $4}' > temp_lj_intra.txt"
    )
    if os.path.exists("temp_lj_intra.txt"):
        try:
            lj_intra_one_seed = np.genfromtxt("temp_lj_intra.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp_lj_intra.txt")
        return np.mean(lj_intra_one_seed[0])
    else:
        return None


def get_tailc():
    """Get tail corection energy for the system."""
    os.system(
        "grep 'Tail correction' run.melt | awk '{print $3}' > temp_tailc.txt"
    )
    if os.path.exists("temp_tailc.txt"):
        try:
            tailc_one_seed = np.genfromtxt("temp_tailc.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp_tailc.txt")
        return np.mean(tailc_one_seed[0])
    else:
        return None


def get_coul():
    """Get coulombic energy for the system."""
    os.system("grep 'coulombic energy' run.melt | awk '{print $3}' > temp.txt")
    if os.path.exists("temp.txt"):
        try:
            one_seed = np.genfromtxt("temp.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp.txt")
        return np.mean(one_seed[0])
    else:
        return None


def get_bond():
    """Get bond energy for the system."""
    os.system("grep 'bond vibration' run.melt | awk '{print $3}' > temp.txt")
    if os.path.exists("temp.txt"):
        try:
            one_seed = np.genfromtxt("temp.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp.txt")
        return np.mean(one_seed[0])
    else:
        return None


def get_angle():
    """Get angle bending energy for the system."""
    os.system("grep 'bond bending' run.melt | awk '{print $3}' > temp.txt")
    if os.path.exists("temp.txt"):
        try:
            one_seed = np.genfromtxt("temp.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp.txt")
        return np.mean(one_seed[0])
    else:
        return None


def get_torsion():
    """Get torsional energy for the system."""
    os.system("grep 'torsional' run.melt | awk '{print $2}' > temp.txt")
    if os.path.exists("temp.txt"):
        try:
            one_seed = np.genfromtxt("temp.txt")
        except Exception as e:
            raise e
        finally:
            # now delete the temp file
            os.remove("temp.txt")
        return np.mean(one_seed[-1])
    else:
        return None


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


if __name__ == "__main__":
    main()
