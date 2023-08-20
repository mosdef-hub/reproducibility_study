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
from matplotlib.ticker import MaxNLocator
from scipy import interpolate, stats
from scipy.stats import norm


def main():
    """Read run files from all the independent seeds and get the final density and RDFs."""
    data_path = "rdf_comparison_data"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    # delete the folder manually
    base_dir = os.getcwd()
    os.chdir(data_path)
    file_flex = {}
    file_rig = {}

    temps = [280, 300, 320]

    for temp in temps:
        file_flex[temp] = (
            base_dir
            + "/rdf_analysis_data/ethanolAA_NPT-flexOH_{:.1f}K_101.325kPa_cutoff_hard_lrc_energy_pressure/avg_rdf.txt".format(
                temp
            )
        )
        file_rig[temp] = (
            base_dir
            + "/rdf_analysis_data/ethanolAA_NPT-fixOH_{:.1f}K_101.325kPa_cutoff_hard_lrc_energy_pressure/avg_rdf.txt".format(
                temp
            )
        )

    rdf_flex = {}
    rdf_rig = {}

    for temp in temps:
        rdf_flex[temp] = np.genfromtxt(file_flex[temp])
        rdf_rig[temp] = np.genfromtxt(file_rig[temp])
        # plot settings
    ms = 8  # markersize
    xtickfs = 11  # xtickfontsize
    xlabelfs = 14  # xlabelfontsize
    ylabelfs = 14  # ylabelfontsize
    ytickfs = 11  # ytickfontsize
    titlefs = 14  # title size
    legendfs = 12
    error_bar_capsize = 3

    sns.set_style("whitegrid")
    fig, axs = plt.subplots(3, 1, sharex=True, sharey=False, figsize=(6, 9))
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    # fig.suptitle('O-O RDF comparison', fontsize = titlefs)

    ax = {}
    ax[0], ax[1], ax[2] = axs[0], axs[1], axs[2]

    for i, temp in enumerate(temps):
        ax[i].plot(
            rdf_flex[temp][:, 0], rdf_flex[temp][:, 1], label="MCCCS-MN-flex"
        )
        ax[i].plot(
            rdf_rig[temp][:, 0], rdf_rig[temp][:, 1], "--", label="MCCCS-MN-fix"
        )

        ax[i].grid(alpha=0.25)
        ax[i].title.set_text(r"$T = {}$".format(temp) + " K")
        ax[i].set_ylim([0, 5])
        ax[i].set_xlim([2, 6])
        ax[i].yaxis.set_major_locator(plt.MaxNLocator(5))
        ax[i].xaxis.set_major_locator(plt.MaxNLocator(4))

    ax[0].legend(frameon=True)
    ax[0].set_ylabel(r"$g(r)$", fontsize=ylabelfs)
    ax[1].set_ylabel(r"$g(r)$", fontsize=ylabelfs)
    ax[2].set_ylabel(r"$g(r)$", fontsize=ylabelfs)
    ax[2].set_xlabel(r"$r$" + r" ($\mathrm{\AA}$)", fontsize=xlabelfs)

    plt.tight_layout()
    plt.savefig("O-O_rdf_comparison.pdf", dpi=500)

    plt.close()

    file_flex = {}
    file_rig = {}

    temps = [280, 300, 320]

    for temp in temps:
        file_flex[temp] = (
            base_dir
            + "/rdf_analysis_data/ethanolAA_NPT-flexOH_{:.1f}K_101.325kPa_cutoff_hard_lrc_energy_pressure/avg_cdf.txt".format(
                temp
            )
        )
        file_rig[temp] = (
            base_dir
            + "/rdf_analysis_data/ethanolAA_NPT-fixOH_{:.1f}K_101.325kPa_cutoff_hard_lrc_energy_pressure/avg_cdf.txt".format(
                temp
            )
        )

    cdf_flex = {}
    cdf_rig = {}

    for temp in temps:
        cdf_flex[temp] = np.genfromtxt(file_flex[temp])
        cdf_rig[temp] = np.genfromtxt(file_rig[temp])
        # plot settings
    ms = 8  # markersize
    xtickfs = 11  # xtickfontsize
    xlabelfs = 14  # xlabelfontsize
    ylabelfs = 14  # ylabelfontsize
    ytickfs = 11  # ytickfontsize
    titlefs = 14  # title size
    legendfs = 12
    error_bar_capsize = 3

    sns.set_style("whitegrid")
    fig, axs = plt.subplots(3, 1, sharex=True, sharey=False, figsize=(6, 9))
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    # fig.suptitle('O-O CDF comparison', fontsize = titlefs)

    ax = {}
    ax[0], ax[1], ax[2] = axs[0], axs[1], axs[2]

    for i, temp in enumerate(temps):
        ax[i].plot(
            rdf_flex[temp][:, 0], cdf_flex[temp][:, 1], label="MCCCS-MN-flex"
        )
        ax[i].plot(
            rdf_rig[temp][:, 0],
            cdf_rig[temp][:, 1],
            "--",
            label="MCCCS-MN-fix",
        )

        f = interpolate.interp1d(rdf_flex[temp][:, 0], cdf_flex[temp][:, 1])
        g = interpolate.interp1d(rdf_rig[temp][:, 0], cdf_rig[temp][:, 1])

        diff = f(3.25) - g(3.25)

        ax[i].grid(alpha=0.25)
        ax[i].title.set_text(r"$T = {}$".format(temp) + " K")
        ax[i].set_ylim([0, 2.4])
        ax[i].set_xlim([2, 4])
        ax[i].yaxis.set_major_locator(plt.MaxNLocator(5))
        ax[i].xaxis.set_major_locator(plt.MaxNLocator(4))
        ax[i].text(
            0.8,
            0.1,
            r"$\Delta n_{r=3.25\;\mathrm{\AA}}=$" + "{:.2f}".format(diff),
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax[i].transAxes,
        )

    ax[0].legend(frameon=True)
    ax[0].set_ylabel(r"$n(r)$", fontsize=ylabelfs)
    ax[1].set_ylabel(r"$n(r)$", fontsize=ylabelfs)
    ax[2].set_ylabel(r"$n(r)$", fontsize=ylabelfs)
    ax[2].set_xlabel(r"$r$" + r" ($\mathrm{\AA}$)", fontsize=xlabelfs)

    plt.tight_layout()
    plt.savefig("O-O_cdf_comparison.pdf", dpi=500)

    plt.close()


if __name__ == "__main__":
    main()
