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
from scipy.signal import argrelextrema


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
        if "NPT" in ensemble and molecule == "methaneUA" and engine == "mcccs":
            rdf_list = []
            cdf_list = []
            for job in group:
                os.chdir(job.ws)
                print(job)
                print(os.listdir())
                rdf_pair = np.genfromtxt("A-A_rdf.txt")
                cdf_pair = np.genfromtxt("A-A_cdf.txt")
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

            # printing the local maxima and minima
            maxima_indices = argrelextrema(avg_rdf[:, 1], np.greater)
            maxima_rs = avg_rdf[:, 0][maxima_indices]

            minima_indices = argrelextrema(avg_rdf[:, 1], np.less)
            minima_rs = avg_rdf[:, 0][minima_indices]

            ones_indices = np.where(
                np.isclose(avg_rdf[:, 1], 1, atol=1e-4) == True
            )[0]
            ones_rs = avg_rdf[:, 0][ones_indices]

            maxima_rdfs = avg_rdf[:, 1][maxima_indices]
            minima_rdfs = avg_rdf[:, 1][minima_indices]
            ones_rdfs = avg_rdf[:, 1][ones_indices]

            maxima_rs = [11.488, 15.0995]
            maxima_rdfs = [1.05, 1.01582]
            minima_rs = [13.2365, 17.0315]
            minima_rdfs = [0.971853, 0.990829]

            ones_rs = [10.66, 12.4545, 14.2715, 16.0425, 17.905]
            ones_rdfs = [1] * 5

            fig, axs = plt.subplots(
                1, 1, sharex=False, sharey=False, figsize=(4.5, 4.5)
            )
            ax1 = axs
            plt.plot(avg_rdf[:, 0], avg_rdf[:, 1])
            plt.grid(alpha=0.25)
            plt.xlabel(r"$r$" + r" ($\mathrm{\AA}$)")
            plt.ylabel(r"$g(r)$")
            plt.title("RDF")
            plt.xlim([3, 18])
            plt.ylim([0, 2.5])

            # Maxima _minima scatters

            plt.scatter(maxima_rs, maxima_rdfs, color="b")

            plt.scatter(minima_rs, minima_rdfs, color="r")

            plt.scatter(ones_rs, ones_rdfs, color="g")

            # adjusting ticks

            for ax in [ax1]:
                ax.xaxis.set_minor_locator(
                    matplotlib.ticker.AutoMinorLocator(2)
                )
            ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))

            plt.tight_layout()

            plt.savefig("A-A_rdf.pdf")
            plt.close()

            fig, axs = plt.subplots(
                1, 1, sharex=False, sharey=False, figsize=(4.5, 4.5)
            )
            ax1 = axs
            plt.plot(avg_cdf[:, 0], avg_cdf[:, 1])
            plt.grid(alpha=0.25)
            plt.xlabel(r"$r$" + r" ($\mathrm{\AA}$)")
            plt.ylabel(r"$n(r)$")
            plt.title("A-A CDF")
            # adjusting ticks

            for ax in [ax1]:
                ax.xaxis.set_minor_locator(
                    matplotlib.ticker.AutoMinorLocator(2)
                )
            ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))
            plt.tight_layout()
            plt.savefig("A-A_cdf.pdf")
            plt.close()

        os.chdir("..")


if __name__ == "__main__":
    main()
