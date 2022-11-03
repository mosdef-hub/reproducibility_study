"""Script for analysing the results from the 16 different seeds from the 11 systems."""
# It also parses the gsd format trajectory stored in each output analysis folder (obtained by executing conv_traj.py before this script) to get the RDFs."""
import os
import shutil
from glob import glob

import freud
import matplotlib
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import seaborn as sns
import signac
from scipy import stats


def main():
    """Read run files from all the independent seeds and get the final density and RDFs."""
    data_path = "hbond_analysis_data"
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
        if "NPT" in ensemble and molecule == "ethanolAA":
            map_output_list = []
            r_list = []
            theta_list = []
            n_hbond_list = []

            for job in group:
                print("job = ", job)
                os.chdir(job.ws)

                map_output = np.loadtxt("map_output.csv", delimiter=",")
                r = np.loadtxt("r.csv", delimiter=",") * 10
                theta = np.loadtxt("theta.csv", delimiter=",")
                n_hbond = np.genfromtxt("n_hbond.txt")
                map_output_list.append(map_output)
                r_list.append(r)
                theta_list.append(theta)
                print(n_hbond)
                n_hbond_list.append(n_hbond)

            os.chdir(base_dir)

            avg_map_output = np.mean(map_output_list, axis=0)
            avg_r = np.mean(r_list, axis=0)
            avg_theta = np.mean(theta_list, axis=0)
            avg_n_hbond = np.mean(n_hbond_list)
            sem = np.std(n_hbond_list) / (len(n_hbond_list)) ** 0.5
            # plot settings
            ms = 8  # markersize
            xtickfs = 11  # xtickfontsize
            xlabelfs = 14  # xlabelfontsize
            ylabelfs = 14  # ylabelfontsize
            ytickfs = 11  # ytickfontsize
            titlefs = 14  # title size
            legendfs = 12
            error_bar_capsize = 3

            plt.figure()
            sns.set_style("whitegrid")
            cmap = plt.get_cmap("jet")
            plt.figure(figsize=(7, 5))
            levels = np.linspace(1, 51, 11)

            cs = plt.contourf(
                avg_r, avg_theta, avg_map_output, levels=levels, cmap=cmap
            )
            plt.xlabel(r"$r$" + r" ($\mathrm{\AA}$)", fontsize=xlabelfs)
            plt.ylabel("\u03B8 (degrees)", fontsize=ylabelfs)
            plt.xlim([2.2, 3.5])
            plt.ylim([125, 180])
            cbar = plt.colorbar(cs)
            cs.axes.tick_params(labelsize=14)

            tick_font_size = 12
            cbar.ax.tick_params(labelsize=tick_font_size)
            ellipse = matplotlib.patches.Ellipse(
                (2.75, 180),
                0.5 * 2,
                50 * 2,
                ec="magenta",
                facecolor="none",
                linewidth=3,
            )
            plt.gca().add_patch(ellipse)
            plt.grid(alpha=0.25)

            plt.tight_layout()
            plt.text(
                3.2,
                130,
                r"$n_\mathrm{HB}=$"
                + "{:.1f}".format(avg_n_hbond)
                + r"$\pm$"
                + "{:.1f}".format(sem),
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=14,
            )

            plt.savefig("hbond.png", dpi=500)
            np.savetxt("nhbond.txt", np.vstack((avg_n_hbond, sem)).T)
            plt.close()

        os.chdir("..")


if __name__ == "__main__":
    main()
