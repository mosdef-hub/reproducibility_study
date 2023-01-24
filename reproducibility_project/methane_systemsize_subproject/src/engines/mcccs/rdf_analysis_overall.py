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
import matplotlib
matplotlib.use("pdf")
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import StrMethodFormatter, NullFormatter, MultipleLocator


from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times']})
#rc('text', usetex=True)


def main():
    """Read run files from all the independent seeds and get the final density and RDFs."""
    data_path = "rdf_analysis_data"
    if os.path.exists(data_path):
        shutil.rmtree(data_path)
    os.makedirs(data_path)
    rdfs = {}
    cdfs = {}
    ensemble2N = {}
    ensemble2N["NPT-exsmall"] = 450
    ensemble2N["NPT-small"] = 600
    ensemble2N["NPT-medium"] = 900
    ensemble2N["NPT-large"] = 1800
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
        ensemble
    ), group in project.groupby(
        (
            "molecule",
            "ensemble",
            "temperature",
            "pressure",
            "cutoff_style",
            "long_range_correction",
            "ensemble",
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
            ensemble,
        )
        if not os.path.isdir(
            "{}_{}_{}K_{}kPa_cutoff_{}_lrc_{}".format(
                molecule,
                ensemble,
                temperature,
                pressure,
                cutoff_style,
                str(long_range_correction),
                ensemble,
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
                    ensemble,
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
                ensemble,
            )
        )

        base_dir = os.getcwd()
        if "NPT" in ensemble and molecule == "methaneUA" and ensemble != "NPT-exlarge":
            rdf_list = []
            cdf_list = []
            for job in group:
                os.chdir(job.ws)
                rdf_pair = np.genfromtxt("C-C_rdf.txt")
                cdf_pair = np.genfromtxt("C-C_cdf.txt")
                rdf_list.append(rdf_pair)
                cdf_list.append(cdf_pair)

            os.chdir(base_dir)
            avg_rdf = np.mean(rdf_list, axis=0)
            avg_rdf[:, 0] = avg_rdf[:, 0] * 10
            np.savetxt("avg_rdf.txt", avg_rdf)

            rdfs[ensemble] = avg_rdf

            avg_cdf = np.mean(cdf_list, axis=0)
            avg_cdf[:, 0] = avg_cdf[:, 0] * 10
            np.savetxt("avg_cdf.txt", avg_cdf)

            cdfs[ensemble] = avg_cdf
            plt.plot(avg_rdf[:, 0], avg_rdf[:, 1])
            plt.grid(alpha=0.25)
            plt.xlabel(r"$r$" + r" ($\mathrm{\AA}$)")
            plt.ylabel(r"$g(r)$")
            plt.title("C-C RDF")
            plt.savefig("C-C_rdf.pdf")
            plt.close()

            plt.plot(avg_cdf[:, 0], avg_cdf[:, 1])
            plt.grid(alpha=0.25)
            plt.xlabel(r"$r$" + r" ($\mathrm{\AA}$)")
            plt.ylabel(r"$n(r)$")
            plt.title("C-C CDF")
            plt.savefig("C-C_cdf.pdf")
            plt.close()
        os.chdir("..")

    fig, axs = plt.subplots(1, 1, sharex=False, sharey=False,figsize=(3,3))
    ax1 = axs
    for key in ["NPT-exsmall", "NPT-small", "NPT-medium", "NPT-large"]:
        ax1.plot(rdfs[key][:,0], rdfs[key][:,1], label = "$N={}$".format(ensemble2N[key]))
    ax1.set_xlabel(r"$r$" + r" ($\mathrm{\AA}$)")
    ax1.set_ylabel(r"$g(r)$")
    ax1.set_xlim([2,12])
    ax1.set_ylim([0.6,1.2])

    plt.title("C-C RDF")
    plt.tight_layout()
    plt.legend()
    #adjusting ticks

    for ax  in [ax1]:
        ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))


    plt.savefig("C-C_rdf_comparison.pdf")


if __name__ == "__main__":
    main()
