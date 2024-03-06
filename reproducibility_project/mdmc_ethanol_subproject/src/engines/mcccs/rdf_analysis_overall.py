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
import signac
from matplotlib import rc
from matplotlib.ticker import (
    FormatStrFormatter,
    MaxNLocator,
    MultipleLocator,
    NullFormatter,
    ScalarFormatter,
    StrMethodFormatter,
)
from scipy import stats


def main():
    """Read run files from all the independent seeds and get the final density and RDFs."""
    data_path = "rdf_analysis_data"
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
            rdf_list = []
            cdf_list = []
            for job in group:
                os.chdir(job.ws)
                rdf_pair = np.genfromtxt("O-O_rdf.txt")
                cdf_pair = np.genfromtxt("O-O_cdf.txt")
                rdf_list.append(rdf_pair)
                cdf_list.append(cdf_pair)

            os.chdir(base_dir)
            avg_rdf = np.mean(rdf_list, axis=0)
            avg_rdf[:, 0] = avg_rdf[:, 0] * 10
            # print(avg_rdf)
            np.savetxt("avg_rdf.txt", avg_rdf)

            avg_cdf = np.mean(cdf_list, axis=0)
            avg_cdf[:, 0] = avg_cdf[:, 0] * 10
            np.savetxt("avg_cdf.txt", avg_cdf)

            fig, axs = plt.subplots(
                1, 1, sharex=False, sharey=False, figsize=(4.5, 4.5)
            )
            ax1 = axs
            plt.plot(avg_rdf[:, 0], avg_rdf[:, 1])
            plt.grid(alpha=0.25)
            plt.xlabel(r"$r$" + r" ($\mathrm{\AA}$)")
            plt.ylabel(r"$g(r)$")
            plt.title("O-O RDF")
            # adjusting ticks

            for ax in [ax1]:
                ax.xaxis.set_minor_locator(
                    matplotlib.ticker.AutoMinorLocator(2)
                )
            ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))

            plt.tight_layout()
            plt.savefig("O-O_rdf.pdf")
            plt.close()

            fig, axs = plt.subplots(
                1, 1, sharex=False, sharey=False, figsize=(4.5, 4.5)
            )
            ax1 = axs
            plt.plot(avg_cdf[:, 0], avg_cdf[:, 1])
            plt.grid(alpha=0.25)
            plt.xlabel(r"$r$" + r" ($\mathrm{\AA}$)")
            plt.ylabel(r"$n(r)$")
            plt.title("O-O CDF")
            # adjusting ticks

            for ax in [ax1]:
                ax.xaxis.set_minor_locator(
                    matplotlib.ticker.AutoMinorLocator(2)
                )
            ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))
            plt.tight_layout()
            plt.savefig("O-O_cdf.pdf")
            plt.close()

        os.chdir("..")


if __name__ == "__main__":
    main()
