"""Script for analysing the results from the 16 different seeds from the 11 systems."""
# It also parses the gsd format trajectory stored in each output analysis folder (obtained by executing conv_traj.py before this script) to get the RDFs."""
import os
import shutil
from glob import glob

import freud
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import signac
from scipy import stats
from scipy.stats import gaussian_kde, norm

# plot settings
ms = 8  # markersize
xtickfs = 11  # xtickfontsize
xlabelfs = 14  # xlabelfontsize
ylabelfs = 14  # ylabelfontsize
ytickfs = 11  # ytickfontsize
titlefs = 14  # title size
legendfs = 9
alpha = 0.2


def main():
    """Read run files from all the independent seeds and get the final density and RDFs."""
    data_path = "bl_analysis_data"
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
        if "NPT" in ensemble and molecule == "ethanolAA":

            bl_list = {}
            mean_bl_list = {}
            mean_bl_list["C-C"] = 0.1529
            mean_bl_list["C-O"] = 0.141
            mean_bl_list["O-H"] = 0.0945
            mean_bl_list["C-H"] = 0.109
            bl_list["C-C"] = []
            bl_list["C-O"] = []
            bl_list["O-H"] = []
            bl_list["C-H"] = []

            for job in group:
                print(job)
                os.chdir(job.ws)
                atom_pairs = [["C", "C"], ["C", "O"], ["O", "H"], ["C", "H"]]

                for pair in atom_pairs:
                    atom1 = pair[0]
                    atom2 = pair[1]
                    bls = np.genfromtxt("bl_{}-{}.txt".format(atom1, atom2))

                    bl_list["{}-{}".format(atom1, atom2)].extend(bls)
            os.chdir(base_dir)
            for pair in atom_pairs:
                atom1 = pair[0]
                atom2 = pair[1]
                sns.set_style("whitegrid")
                plt.figure(figsize=(4.5, 4.5))
                numbers = bl_list["{}-{}".format(atom1, atom2)]
                kde = gaussian_kde(numbers)
                nbins = 100

                # Evaluate the KDE on a grid of points
                x_grid = np.linspace(min(numbers), max(numbers), nbins)
                kde_values = kde.evaluate(x_grid)
                # Plot the KDE as a smooth curve
                fig, ax = plt.subplots(figsize=(4.5, 4.5))
                n, bins, patches = ax.hist(
                    x_grid,
                    bins=nbins,
                    weights=kde_values,
                    color="#1f77b4",
                    alpha=0.7,
                )

                # Stack the freq and bin centers vertically
                bin_centers = (bins[:-1] + bins[1:]) / 2
                data = np.column_stack((bin_centers, n))

                plt.axvline(
                    mean_bl_list["{}-{}".format(atom1, atom2)] * 10,
                    ls="--",
                    color="r",
                    label="FF eq.\nbond length",
                )
                plt.legend(
                    frameon=True,
                    loc="upper right",
                    ncol=1,
                    fontsize=legendfs,
                    labelspacing=0.05,
                )

                plt.xlabel(
                    "Bond length" + r" ($\mathrm{\AA}$)", fontsize=xlabelfs
                )
                plt.ylabel("Probability", fontsize=ylabelfs)
                plt.grid(alpha=0.25)
                plt.xticks(fontsize=xtickfs)
                plt.yticks(fontsize=ytickfs)
                plt.title("{}-{}".format(atom1, atom2), fontsize=titlefs)
                plt.tight_layout()

                plt.savefig("{}-{}.pdf".format(atom1, atom2))
                np.savetxt("{}_{}_hist.txt".format(atom1, atom2), data)
                plt.close()

            fig, axs = plt.subplots(
                2, 2, sharex=False, sharey=False, figsize=(9, 9)
            )
            fig.subplots_adjust(hspace=0.3)
            fig.subplots_adjust(wspace=0.3)
            # fig.suptitle(
            #    "Bond length distribution ethanol-AA {} K".format(temperature),
            #    fontsize=titlefs,
            # )

            ax = {}
            ax[0], ax[1], ax[2], ax[3] = (
                axs[0, 0],
                axs[0, 1],
                axs[1, 0],
                axs[1, 1],
            )

            for i, pair in enumerate(atom_pairs):
                atom1 = pair[0]
                atom2 = pair[1]

                numbers = bl_list["{}-{}".format(atom1, atom2)]
                kde = gaussian_kde(numbers)
                nbins = 100

                # Evaluate the KDE on a grid of points
                x_grid = np.linspace(min(numbers), max(numbers), nbins)
                kde_values = kde.evaluate(x_grid)

                a = ax[i].hist(
                    x_grid,
                    bins=nbins,
                    weights=kde_values,
                    color="#1f77b4",
                    ec="#1f77b4",
                    alpha=0.5,
                )
                sns.set_style("whitegrid")
                hist = a[0]
                bins = a[1]
                center = (bins[:-1] + bins[1:]) / 2
                mean, std = norm.fit(bl_list["{}-{}".format(atom1, atom2)])
                # print("The mean is {} and std is {}".format(mean, std))
                # print(center)
                if i == 0:
                    ax[i].axvline(
                        mean_bl_list["{}-{}".format(atom1, atom2)] * 10,
                        ls="--",
                        color="r",
                        label="FF eq.\nbond length",
                    )
                else:
                    ax[i].axvline(
                        mean_bl_list["{}-{}".format(atom1, atom2)] * 10,
                        ls="--",
                        color="r",
                    )

                if i == 2:
                    y = norm.pdf(center, mean, std)
                    ax[i].plot(
                        center,
                        y,
                        c="#ff7f0e",
                        alpha=0.7,
                        label="Gaussian fit\n"
                        + r"$r_0$="
                        + "{:.3f} ".format(mean)
                        + r"$\mathrm{\AA}$",
                    )
                    diff = mean - 0.945

                ax[i].grid(alpha=0.25)
                ax[i].title.set_text("{}-{}".format(atom1, atom2))
            print("mean is {}".format(mean))
            ax[2].text(
                0.2,
                0.7,
                r"$\Delta r_0={:.3f}$".format(diff) + " " + r"$\mathrm{\AA}$",
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax[2].transAxes,
            )
            ax[0].legend(frameon=True)
            ax[2].legend(frameon=True)
            ax[0].set_ylabel("Probability density", fontsize=ylabelfs)
            ax[2].set_ylabel("Probability density", fontsize=ylabelfs)
            ax[2].set_xlabel(
                "Bond length" + r" ($\mathrm{\AA}$)", fontsize=xlabelfs
            )
            ax[3].set_xlabel(
                "Bond length" + r" ($\mathrm{\AA}$)", fontsize=xlabelfs
            )
            ax[0].set_xlim([1.4, 1.65])
            ax[1].set_xlim([1.3, 1.55])
            ax[2].set_xlim([0.85, 1.05])
            ax[3].set_xlim([0.95, 1.2])
            plt.tight_layout()
            plt.savefig("bl_T={}.pdf".format(temperature))
            plt.close()

        os.chdir("..")


if __name__ == "__main__":
    main()
