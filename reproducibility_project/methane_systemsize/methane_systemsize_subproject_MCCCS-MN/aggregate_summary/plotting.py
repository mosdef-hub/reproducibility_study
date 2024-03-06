"""Plotting code for methane system size subproject."""

import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.stats
import seaborn as sns
from matplotlib import container

from reproducibility_project.src.utils.plotting import (
    colors,
    figsize,
    mean_confidence_interval,
    pretty_names,
    symbols,
    wrap_labels,
)


def data_cleaning(molecule_set, summary_df):
    """Make sure properties of all molecule in the provided molecule_set has been calculated."""
    for molecule in molecule_set:
        try:
            mol_group = summary_df.groupby(group_key)
            mol_df = mol_group.get_group(molecule)
        except KeyError:
            print(f"skipping: {molecule}, no data available.")
            continue
        mol_df["statepoint"] = (
            summary_df["temperature"].map(str)
            + "K, "
            + summary_df["pressure"].map(str)
            + "kPa, "
            + summary_df["ensemble"]
        )
        sp_set = set(mol_df["statepoint"])
        goods = list()
        bads = list()
        for engine in mol_df["engine"].unique():
            found_sp = mol_df[mol_df["engine"] == engine].index.tolist()
            if len(found_sp) == len(set(mol_df["statepoint"])):
                goods.append((engine, found_sp))
            else:
                mol_df.drop(axis=0, labels=found_sp, inplace=True)
                bads.append((engine, found_sp))
        print(molecule)
        if goods:
            print(f"\tPassed: {[good[0] for good in goods]}")
        if bads:
            print(f"\tFailed: {[bad[0] for bad in bads]}")


def create_deviation_with_sem_plots(
    molecules,
    engines,
    figsize,
    summary_df,
    group_key,
    colors=colors,
    pretty_names=pretty_names,
    symbols=symbols,
):
    """Create deviation plots of provided molecules and engines lists."""
    for molecule in molecules:
        fig, ax = plt.subplots(figsize=figsize)
        try:
            mol_group = summary_df.groupby(group_key)
            mol_df = mol_group.get_group(molecule)
        except KeyError:
            print(f"skipping: {molecule}, no data available.")
            continue

        mol_df["statepoint"] = (
            summary_df["liqN"].map(str)
            + "\n"
            + summary_df["temperature"].map(str)
            + "K, "
            + summary_df["pressure"].map(str)
            + "kPa"
        )
        mol_df.sort_values(by="liqN", inplace=True)

        statepoints = list(set(mol_df["statepoint"]))
        statepoints_overall_means = list()
        xticks = list()
        sps_positions = list()
        for n_statepoint, statepoint in enumerate(statepoints):
            densities = list()
            stds = list()
            sems = list()
            Ns = list()
            for engine in engines:
                shortlisted = mol_df[
                    (mol_df["engine"] == engine)
                    & (mol_df["statepoint"] == statepoint)
                ]
                densities.append(shortlisted["density-avg"])
                stds.append(shortlisted["density-std"])
                sems.append(shortlisted["density-sem"])
                Ns.append(shortlisted["liqN"])
            overall_mean = np.mean(densities)
            statepoints_overall_means.append(overall_mean)

            confidence_interval = mean_confidence_interval(
                densities, np.array(sems), confidence=0.95
            )
            percentage_delta_density = (
                (densities - overall_mean) * 100 / overall_mean
            )
            percentage_sem_density = 100 * np.array(sems) / overall_mean
            percentage_ci_density = (
                100 * np.array(confidence_interval) / overall_mean
            )

            ind = list()
            sp_position = list()
            for i, engine in enumerate(engines):
                ind.append(i)
                # sp_position.append(1*n_statepoint + 0.1*i)
                sp_position.append(Ns[i])
                ax.errorbar(
                    # 1*n_statepoint + 0.1*i,
                    Ns[i],
                    percentage_delta_density[i],
                    marker=symbols[engine],
                    yerr=percentage_ci_density[i],
                    color=colors[engine],
                    ls="",
                    label=engine,
                )

            sps_positions.append(np.mean(sp_position))
            xticks.append(
                statepoint
                + r"$\rho_{\mathrm{ave}}$"
                + r"$ = {:.4f}$".format(overall_mean)
            )
        ax.set_xlabel(r"N")
        ax.set_ylabel(r"$\frac{100\times\Delta\rho}{\rho}$")
        ax.tick_params(axis="y")
        # plt.title(f"{molecule}")
        props = dict(boxstyle="round", facecolor="none", alpha=1, ec="grey")
        # string+='{:.5f}'.format(overall_mean)

        ax.set_xticks([pos for pos in sps_positions])
        # ax.set_xticklabels([
        #     f'{sp}\n' +r'$\rho_{\mathrm{ave}}$'+'$ = {:.4f}$'.format(sp_omean)
        #     for sp_omean, sp in zip(statepoints_overall_means, statepoints)
        #     ]
        # )

        # Handle ylim
        low, high = ax.get_ylim()
        bound = max(abs(low), abs(high))
        ax.set_ylim(-bound * 1.1, bound * 1.1)

        # get handles
        # wrap_labels(ax, 10)
        handles, labels = ax.get_legend_handles_labels()

        # Sorting handles and labels:
        labels_handles = dict()
        for handle, label in zip(handles, labels):
            labels_handles[label] = handle
        sorted_handles = [labels_handles[engine] for engine in engines]
        sorted_labels = [pretty_names[engine] for engine in engines]
        plt.legend(
            sorted_handles,
            sorted_labels,
            facecolor="white",
            loc="best",
            ncol=2,
            prop={"size": 12},
        )
        # wrap_labels(ax, 10)
        plt.tight_layout()
        plt.grid(alpha=0.0, axis="x")
        plt.savefig(f"figures/{molecule}_deviation_with_sem.pdf", dpi=500)

        plt.close()


def create_density_with_std_plots(
    molecules,
    engines,
    figsize,
    summary_df,
    group_key,
    colors=colors,
    pretty_names=pretty_names,
    symbols=symbols,
):
    """Create density with std plots of provided molecules and engines lists."""
    for molecule in molecules:
        fig, ax = plt.subplots(figsize=figsize)
        try:
            mol_group = summary_df.groupby(group_key)
            mol_df = mol_group.get_group(molecule)
        except KeyError:
            print(f"skipping: {molecule}, no data available.")
            continue

        mol_df["statepoint"] = (
            summary_df["liqN"].map(str)
            + "\n"
            + summary_df["temperature"].map(str)
            + "K, "
            + summary_df["pressure"].map(str)
            + "kPa"
        )
        mol_df.sort_values(by="liqN", inplace=True)

        statepoints = list(set(mol_df["statepoint"]))
        statepoints_overall_means = list()

        xticks = list()
        sps_positions = list()
        for n_statepoint, statepoint in enumerate(statepoints):
            densities = list()
            stds = list()
            sems = list()
            Ns = list()
            for engine in engines:
                shortlisted = mol_df[
                    (mol_df["engine"] == engine)
                    & (mol_df["statepoint"] == statepoint)
                ]
                densities.append(shortlisted["density-avg"])
                stds.append(shortlisted["density-std"])
                sems.append(shortlisted["density-sem"])
                Ns.append(shortlisted["liqN"])
            overall_mean = np.mean(densities)
            statepoints_overall_means.append(overall_mean)

            ind = list()
            sp_position = list()
            for i, engine in enumerate(engines):
                ind.append(i)
                # sp_position.append(1*n_statepoint + 0.1*i)
                sp_position.append(Ns[i])
                ax.errorbar(
                    # 1*n_statepoint + 0.1*i,
                    Ns[i],
                    densities[i],
                    marker=symbols[engine],
                    yerr=stds[i],
                    color=colors[engine],
                    ls="",
                    label=engine,
                )

            sps_positions.append(np.mean(sp_position))
            xticks.append(
                statepoint
                + r"$\rho_{\mathrm{ave}}$"
                + r"$ = {:.4f}$".format(overall_mean)
            )
        ax.set_xlabel("N")
        ax.set_ylabel(r"$\frac{100\times\Delta\rho}{\rho}$")
        ax.tick_params(axis="y")
        # plt.title(f"{molecule}")
        props = dict(boxstyle="round", facecolor="none", alpha=1, ec="grey")
        # string+='{:.5f}'.format(overall_mean)

        ax.set_xticks([pos for pos in sps_positions])

        handles, labels = ax.get_legend_handles_labels()

        # Sorting handles and labels:
        labels_handles = dict()
        for handle, label in zip(handles, labels):
            labels_handles[label] = handle
        sorted_handles = [labels_handles[engine] for engine in engines]
        sorted_labels = [pretty_names[engine] for engine in engines]
        plt.legend(
            sorted_handles,
            sorted_labels,
            facecolor="white",
            loc="best",
            ncol=2,
            prop={"size": 12},
        )
        plt.tight_layout()
        plt.grid(alpha=0.0, axis="x")
        plt.savefig(f"figures/{molecule}_density_with_std.pdf", dpi=500)

        plt.close()


def create_density_with_sem_plots(
    molecules,
    engines,
    figsize,
    summary_df,
    group_key,
    colors=colors,
    pretty_names=pretty_names,
    symbols=symbols,
):
    """Create density with sem plots of provided molecules and engines lists."""
    for molecule in molecules:
        fig, ax = plt.subplots(figsize=figsize)
        try:
            mol_group = summary_df.groupby(group_key)
            mol_df = mol_group.get_group(molecule)
        except KeyError:
            print(f"skipping: {molecule}, no data available.")
            continue

        mol_df["statepoint"] = (
            summary_df["liqN"].map(str)
            + "\n"
            + summary_df["temperature"].map(str)
            + "K, "
            + summary_df["pressure"].map(str)
            + "kPa"
        )
        mol_df.sort_values(by="liqN", inplace=True)

        statepoints = list(set(mol_df["statepoint"]))
        statepoints_overall_means = list()

        xticks = list()
        sps_positions = list()
        for n_statepoint, statepoint in enumerate(statepoints):
            densities = list()
            stds = list()
            sems = list()
            Ns = list()
            for engine in engines:
                shortlisted = mol_df[
                    (mol_df["engine"] == engine)
                    & (mol_df["statepoint"] == statepoint)
                ]
                densities.append(shortlisted["density-avg"])
                stds.append(shortlisted["density-std"])
                sems.append(shortlisted["density-sem"])
                Ns.append(shortlisted["liqN"])
            overall_mean = np.mean(densities)
            statepoints_overall_means.append(overall_mean)

            ind = list()
            sp_position = list()
            for i, engine in enumerate(engines):
                ind.append(i)
                # sp_position.append(1*n_statepoint + 0.1*i)
                sp_position.append(Ns[i])
                ax.errorbar(
                    # 1*n_statepoint + 0.1*i,
                    Ns[i],
                    densities[i],
                    marker=symbols[engine],
                    yerr=sems[i],
                    color=colors[engine],
                    ls="",
                    label=engine,
                )

            sps_positions.append(np.mean(sp_position))
            xticks.append(
                statepoint
                + r"$\rho_{\mathrm{ave}}$"
                + r"$ = {:.4f}$".format(overall_mean)
            )
        ax.set_xlabel("N")
        ax.set_ylabel(r"$\frac{100\times\Delta\rho}{\rho}$")
        ax.tick_params(axis="y")
        # plt.title(f"{molecule}")
        props = dict(boxstyle="round", facecolor="none", alpha=1, ec="grey")

        ax.set_xticks([pos for pos in sps_positions])

        handles, labels = ax.get_legend_handles_labels()

        # Sorting handles and labels:
        labels_handles = dict()
        for handle, label in zip(handles, labels):
            labels_handles[label] = handle
        sorted_handles = [labels_handles[engine] for engine in engines]
        sorted_labels = [pretty_names[engine] for engine in engines]
        plt.legend(
            sorted_handles,
            sorted_labels,
            facecolor="white",
            loc="best",
            ncol=2,
            prop={"size": 12},
        )
        plt.tight_layout()
        plt.grid(alpha=0.0, axis="x")
        plt.savefig(f"figures/{molecule}_density_with_sem.pdf", dpi=500)

        plt.close()


if __name__ == "__main__":
    import os

    if not os.path.isdir("figures"):
        os.mkdir("figures")

    summary_df = pd.read_csv("aggregate_summary_all.csv", index_col=0)
    group_key = "molecule"

    molecules = ("methaneUA",)
    engines = (
        "lammps-VU",
        "mcccs",
    )
    N_mols = [450, 600, 900, 1800]
    ensembles = ["NPT-exsmall", "NPT-small", "NPT-medium", "NPT-large"]

    mol_ensemble_conversion = dict()
    for mol, ensemble in zip(N_mols, ensembles):
        mol_ensemble_conversion[ensemble] = mol

    summary_df["liqN"] = [
        mol_ensemble_conversion[ensemble] for ensemble in summary_df["ensemble"]
    ]

    plt.rcParams["font.family"] = "helvetica"

    data_cleaning(molecules, summary_df)
    create_density_with_sem_plots(
        molecules=molecules,
        engines=engines,
        figsize=figsize,
        summary_df=summary_df,
        group_key=group_key,
    )
    create_density_with_std_plots(
        molecules=molecules,
        engines=engines,
        figsize=figsize,
        summary_df=summary_df,
        group_key=group_key,
    )
    create_deviation_with_sem_plots(
        molecules=molecules,
        engines=engines,
        figsize=figsize,
        summary_df=summary_df,
        group_key=group_key,
    )
