"""Functionality to generate figure plots to compare MoSDeF to RR relative errors."""

import warnings

warnings.filterwarnings("ignore")
import copy
import os
import re
import shutil
from collections import Counter

# Plot settings and imports
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.use("pdf")
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
import numpy as np
from matplotlib import rc
from matplotlib.ticker import (
    AutoMinorLocator,
    FormatStrFormatter,
    MaxNLocator,
    MultipleLocator,
    NullFormatter,
    ScalarFormatter,
    StrMethodFormatter,
)

rc("font", **{"family": "serif", "serif": ["Times"]})
rc("text", usetex=True)

# plot settings
ms = 8  # markersize
xtickfs = 15  # xtickfontsize
xlabelfs = 18  # xlabelfontsize
ylabelfs = 20  # ylabelfontsize
ytickfs = 14  # ytickfontsize
titlefs = 14  # title size
legendfs = 9
alpha = 0.2


# equation to calculate RE for a pandas series
def _calculate_relative_error(ser):
    dataArray = ser.to_numpy()
    mean = np.nanmean(dataArray)
    return (dataArray - mean) / mean * 100


# 2 Load RR data
def _read_rr_csv(model):
    df = pd.read_csv(f"csvs/{model}_data.csv", index_col=0)
    pressure = float(model.split(" ")[3])
    df.insert(len(df.columns) - 2, "pressure", np.full(len(df.index), pressure))
    return df


# 4 filter RR data
def _filter_rr_data(df, model, remove_MC=True):
    keep_labelsList = ["Program(Group)", "T / K", "ρ / kg m−3", "pressure"]
    new_labelsList = ["engine", "temperature", "density", "pressure"]
    df = df[keep_labelsList]
    df.columns = new_labelsList
    MCList = ["TOWHEE(BS)", "Ind.MC(PB)" "ms2*(KL)"]

    n_rows = len(df.index)
    df.insert(1, "replicate", np.full(n_rows, 0), True)
    df.insert(4, "density-std", np.full(n_rows, 0.0), True)
    df.insert(5, "forcefield", np.full(n_rows, model.split()[0]), True)
    df.insert(1, "molecule", np.full(n_rows, model.split()[1]), True)
    df.insert(0, "associated_work", np.full(n_rows, "RR"), True)
    df.insert(
        1,
        "MCorMD",
        list(map(lambda x: "MC" if x in MCList else "MD", df["engine"])),
    )

    return df


def _load_all_rr_data():
    rr_moleculeList = [
        "OPLS Ethane at 5 MPa",
        "OPLS Ethane at 41 MPa",
        "OPLS Ethane at 70 MPa",
        "TraPPE Ethane at 5 MPa",
        "TraPPE Ethane at 41 MPa",
        "TraPPE Ethane at 70 MPa",
        "OPLSAMBER Ethane at 5 MPa",
        "OPLSAMBER Ethane at 41 MPa",
        "OPLSAMBER Ethane at 70 MPa",
        "OPLS Propane at 5 MPa",
        "OPLS Propane at 41 MPa",
        "OPLS Propane at 70 MPa",
        "TraPPE Propane at 5 MPa",
        "TraPPE Propane at 41 MPa",
        "TraPPE Propane at 70 MPa",
        "OPLSAMBER Propane at 5 MPa",
        "OPLSAMBER Propane at 41 MPa",
        "OPLSAMBER Propane at 70 MPa",
        "OPLS n-Butane at 5 MPa",
        "OPLS n-Butane at 41 MPa",
        "OPLS n-Butane at 70 MPa",
        "TraPPE n-Butane at 5 MPa",
        "TraPPE n-Butane at 41 MPa",
        "TraPPE n-Butane at 70 MPa",
        "OPLSAMBER n-Butane at 5 MPa",
        "OPLSAMBER n-Butane at 41 MPa",
        "OPLSAMBER n-Butane at 70 MPa",
        "OPLS iso-Butane at 5 MPa",
        "OPLS iso-Butane at 41 MPa",
        "OPLS iso-Butane at 70 MPa",
        "TraPPE iso-Butane at 5 MPa",
        "TraPPE iso-Butane at 41 MPa",
        "TraPPE iso-Butane at 70 MPa",
        "OPLSAMBER iso-Butane at 5 MPa",
        "OPLSAMBER iso-Butane at 41 MPa",
        "OPLSAMBER iso-Butane at 70 MPa",
    ]
    hasse_dfList = []
    for molecule in rr_moleculeList:
        hasseDF = _read_rr_csv(molecule)
        hasse_dfList.append(_filter_rr_data(hasseDF, molecule, remove_MC=False))
    return hasse_dfList


def _mask_df(df):
    # mask = ((df["molecule"] == "Ethane") & (df["forcefield"] == "TraPPE") #testing masks to check for outliers
    #       & (df["temperature"] == 298.0) & (df["engine"] == "ms2(KL)")
    #        & (df["pressure"] == 5)
    #       )
    # df = df[~mask]
    mask = (df["molecule"] == "Ethane") & (df["engine"] == "ms2(KL)")
    df = df[~mask]
    mask = (df["molecule"] == "Ethane") & (df["engine"] == "GROMACS(KL)")
    df = df[~mask]
    # mask = ((df["molecule"] == "Ethane") & (df["forcefield"] == "TraPPE")
    #    & (df["temperature"] == 298.0) & (df["engine"] == "GROMACS(KL)")
    #    & (df["pressure"] == 5)
    #   )
    df = df[~mask]
    df = df[
        ~(
            (df["temperature"] == 298.0)
            & (df["forcefield"] == "OPLSAMBER")
            & (df["molecule"] == "Ethane")
        )
    ]
    df = df[~(df["molecule"] == "pentaneUA-flexible_bonds")]
    df = df[~((df["engine"] == "gromacs") & (df["molecule"] == "benzeneUA"))]

    return df


def _transform_ff_names(df):
    df["forcefield"] = df["forcefield"].transform(
        lambda x: "TraPPE" if "ua" in x else x
    )
    df["forcefield"] = df["forcefield"].transform(
        lambda x: x.upper() if x == "oplsaa" else x
    )
    df["forcefield"] = df["forcefield"].transform(
        lambda x: "SPC/E" if x == "spce" else x
    )
    df["ff-mol"] = df[["molecule", "forcefield"]].agg("-".join, axis=1)
    df["ff-mol"] = df["ff-mol"].transform(
        lambda x: (
            re.sub("([a-zA-Z])", lambda y: y.groups()[0].upper(), x, 1)
            if "Butane" not in x
            else x
        )
    )


def _sort_plotting_labels(df, x_label, y_label):
    # Prepare colors
    sortmolDict = {
        "Ethane": 0,
        "Propane": 1,
        "n-Butane": 2,
        "iso-Butane": 3,
        "methaneUA": 4,
        "pentaneUA-constrain_bonds": 5,
        "pentaneUA-flexible_bonds": 9,
        "benzeneUA": 6,
        "waterSPCE": 7,
        "ethanolAA": 8,
    }
    sortffmolDict = {
        "Ethane-TraPPE": 1,
        "Propane-TraPPE": 2,
        "n-Butane-TraPPE": 3,
        "iso-Butane-TraPPE": 4,
        "Ethane-OPLS": 5,
        "Propane-OPLS": 6,
        "n-Butane-OPLS": 7,
        "iso-Butane-OPLS": 8,
        "Ethane-OPLSAMBER": 9,
        "Propane-OPLSAMBER": 10,
        "n-Butane-OPLSAMBER": 11,
        "iso-Butane-OPLSAMBER": 12,
        "MethaneUA-TraPPE": 13,
        "PentaneUA-constrain_bonds-TraPPE": 14,
        "BenzeneUA-TraPPE": 16,
        "WaterSPCE-SPC/E": 17,
        "EthanolAA-OPLSAA": 18,
    }
    if x_label == "molecule":
        df["sortedindex"] = df["molecule"].apply(
            func=lambda x: sortmolDict.get(x)
        )
    elif x_label == "ff-mol":
        df["sortedindex"] = df[x_label].apply(
            func=lambda x: sortffmolDict.get(x)
        )
    df = df.sort_values(by=["sortedindex"], ascending=True).copy()
    color_labels = list(Counter(df[x_label].to_numpy()).keys())
    palette = sns.color_palette("colorblind", n_colors=len(color_labels))
    color_map = dict(zip(color_labels, palette))

    xtick_labelsList = copy.deepcopy(color_labels)
    for i in range(len(color_labels)):
        for filter_label in [
            "AA",
            "SPCE",
            "UA",
            "-constrain_bonds",
            "-flexible_bonds",
        ]:
            if filter_label in color_labels[i]:
                xtick_labelsList[i] = xtick_labelsList[i].replace(
                    filter_label, ""
                )
    return color_labels, xtick_labelsList, color_map


def _plot_params(color_labels):
    """Generate summary plots for the aggregate project."""
    fig, ax = plt.subplots(
        1,
        1,
        sharex=False,
        sharey=False,
        figsize=(len(color_labels) * 0.5 + 3, 4),
    )
    trans = plt.gca().transData
    x_ticks = np.arange(len(color_labels))
    ax.set_xticks(x_ticks)
    ax.set_axisbelow(True)
    ax.grid(visible=True, which="major", axis="y", alpha=0.5, color="gray")
    ax.set_xticks([], minor=True)
    ax.set_yticks([], minor=True)

    return fig, ax


def _plot_scatter(df, fig, ax, color_labels, color_map, x_label, y_label):
    x_ticks = np.arange(len(color_labels))
    for i, molecule in enumerate(color_labels):
        series = df.loc[df[x_label] == molecule][y_label]
        ydata = series[series.notnull()]
        if len(x_ticks) > 1:
            spread = (x_ticks[1] - x_ticks[0]) * 0.25
        else:
            spread = 0.25
        offsets = np.linspace(-1 * spread, spread, len(ydata))
        xdata = np.full(len(ydata), x_ticks[i]) + offsets
        ax.scatter(
            x=xdata,
            y=ydata,
            color=color_map[molecule],
            label=molecule,
            s=np.full(len(ydata.index), 50),
            edgecolors="black",
            marker="o",
            alpha=0.6,
        )
        # ax.axvline(x=6.5, linestyle="--", color="k")


def _plot_violin(df, fig, ax, color_labels, color_map, x_label, y_label):
    x_ticks = np.arange(len(color_labels))
    for i, molecule in enumerate(color_labels):
        series = df.loc[df[x_label] == molecule][y_label]
        ydata = series[series.notnull()]
        if len(x_ticks) > 1:
            spread = (x_ticks[1] - x_ticks[0]) * 0.25
        else:
            spread = 0.25
        print(
            f"Maximum deviation is: {max(abs(ydata)):.2f} for molecule {color_labels[i]}"
        )
        offsets = np.linspace(-1 * spread, spread, len(ydata))
        xdata = np.full(len(ydata), x_ticks[i]) + offsets
        ax.violinplot([ydata], positions=[i], showextrema=True)


def _plot_box(df, fig, ax, color_labels, color_map, x_label, y_label):
    x_ticks = np.arange(len(color_labels))
    for i, molecule in enumerate(color_labels):
        series = df.loc[df[x_label] == molecule][y_label]
        ydata = series[series.notnull()]
        if len(x_ticks) > 1:
            spread = (x_ticks[1] - x_ticks[0]) * 0.25
        else:
            spread = 0.25
        offsets = np.linspace(-1 * spread, spread, len(ydata))
        xdata = np.full(len(ydata), x_ticks[i]) + offsets
        ax.boxplot([ydata], positions=[i])


import textwrap


def _wrap_labels(ax, width, break_long_words=False, rotation=45):
    labels = []
    for label in ax.get_xticklabels():
        text = label.get_text()
        labels.append(
            textwrap.fill(text, width=width, break_long_words=break_long_words)
        )
    if rotation > 0:
        ax.set_xticklabels(
            labels,
            rotation_mode="anchor",
            rotation=rotation,
            ha="right",
            size=xtickfs,
        )
    else:
        ax.set_xticklabels(
            labels,
            rotation_mode="anchor",
            rotation=rotation,
            ha="center",
            size=xtickfs,
        )


def _general_plot_formatting(
    df,
    ax,
    splitby,
    figdir,
    combMCMD,
    xtick_labelsList,
    color_labels,
    average_replicas,
    ymax=None,
    rotation=0,
    wrap_labels=True,
):
    # tick params
    ax.tick_params(axis="y", labelsize=ytickfs)
    ax.tick_params(axis="x", labelsize=xtickfs)
    ax.set_xlim([-0.5, len(color_labels) - 0.5])
    if ymax is None:
        maxy = max(np.abs(df["Relative_Error"])) * 1.05
    else:
        maxy = ymax
    ax.set_ylim([-maxy, maxy])
    ax.yaxis.set_major_locator(plt.MaxNLocator(3))
    x_ticks = np.arange(len(color_labels))
    ax.set_xticks(x_ticks)
    ax.set_axisbelow(True)
    ax.grid(visible=True, which="major", axis="y", alpha=0.5, color="gray")
    # ax.axhline(y = 0, color = 'black', alpha=0.5, linestyle = '--', zorder=-1)
    ax.tick_params(which="both", width=1)
    ax.set_xticks([], minor=True)
    ax.tick_params(which="minor", length=4, axis="y", bottom=True)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    if splitby == "associated_work":
        ax.set_xticklabels(
            xtick_labelsList, ha="center", rotation_mode="anchor", size=15
        )
    elif figdir == "MoSDeF":
        xtick_labelsList = list(
            map(lambda lab: "-\n".join(lab.split("-")), xtick_labelsList)
        )
        ax.set_xticklabels(xtick_labelsList, size=xlabelfs, rotation=rotation)
    else:
        ax.set_xticklabels(
            xtick_labelsList,
            size=xlabelfs,
            rotation_mode="anchor",
            rotation=rotation,
            ha="right",
        )
    if wrap_labels:
        _wrap_labels(ax, width=12, break_long_words=True, rotation=rotation)
    ax.set_ylabel(r"$100\times\Delta\rho / \rho$", size=ylabelfs)


def _save_single_scatter_labels(
    fig, ax, combMCMD, average_replicas, figdir, splitby, plot_type
):
    if combMCMD:
        MCMDStr = "MCandMD"
    else:
        MCMDStr = "MCorMD"

    if average_replicas:
        avgStr = "_avg"
    else:
        avgStr = ""
    fname = (
        f"figures/{figdir}/{figdir}_{splitby}_{MCMDStr}{avgStr}_{plot_type}.pdf"
    )
    print(f"Plotting {fname}")
    fig.tight_layout()
    fig.savefig(fname, dpi=900)


def single_plot_comparisons(
    df, splitby, combMCMD, plot_type, average_replicas=False, figdir="combined"
):
    """Plot deviations split into splitby groupings on xaxis and compared to their statepoint and/or MD/MC means."""
    if average_replicas:
        df = copy.deepcopy(
            df.groupby(
                [
                    "molecule",
                    "forcefield",
                    "engine",
                    "MCorMD",
                    "associated_work",
                    "temperature",
                    "pressure",
                ]
            )["density"]
            .mean()
            .reset_index()
        )

    df = _mask_df(df)

    if combMCMD:
        groupREList = ["molecule", "temperature", "forcefield", "pressure"]
    else:
        groupREList = [
            "MCorMD",
            "molecule",
            "temperature",
            "forcefield",
            "pressure",
        ]
    df["Relative_Error"] = df.groupby(groupREList)["density"].transform(
        _calculate_relative_error
    )

    _transform_ff_names(df)

    x_label = splitby  #'molecule', "forcefield", "associated_work"
    y_label = "Relative_Error"

    color_labels, xtick_labelsList, color_map = _sort_plotting_labels(
        df, x_label, y_label
    )
    fig, ax = _plot_params(color_labels)

    if plot_type == "scatter":
        _plot_scatter(df, fig, ax, color_labels, color_map, x_label, y_label)
    elif plot_type == "violin":
        _plot_violin(df, fig, ax, color_labels, color_map, x_label, y_label)
    else:
        return

    _general_plot_formatting(
        df,
        ax,
        splitby,
        figdir,
        groupMCMD,
        xtick_labelsList,
        color_labels,
        average_replicas,
    )  # formatting for all plots
    _save_single_scatter_labels(
        fig, ax, groupMCMD, average_replicas, figdir, splitby, plot_type
    )
    plt.close()

    return ax


def _plot_stack_params(stack):
    """Initialize subplots with shared x axis."""
    fig = plt.figure(figsize=(9, 6))
    gs = fig.add_gridspec(len(stack), hspace=0)
    stackax = gs.subplots(sharex=True, sharey=False)

    return fig, stackax


def _plot_stack_params_multix(stack, hspace=0.5):
    """Initialize subplots with shared x axis."""
    fig = plt.figure(figsize=(9, 9))
    gs = fig.add_gridspec(len(stack), hspace=hspace)
    stackax = gs.subplots(sharex=False, sharey=False)

    return fig, stackax


def stack_plot_comparisons(df, stack, ymax=None):
    """Stack multiple plots on each other to generate a single figure."""
    # Initialize Figure
    fig, axs = _plot_stack_params(stack)

    # Iterate through stack to generate plot
    for ax, sta in zip(axs, stack):
        if sta["figdir"] in ["MoSDeF", "RR"]:
            stackDF = df.loc[df["associated_work"] == sta["figdir"]].copy()
        else:
            stackDF = df.copy()

        if sta["average_replicas"]:
            stackDF = copy.deepcopy(
                stackDF.groupby(
                    [
                        "molecule",
                        "forcefield",
                        "engine",
                        "MCorMD",
                        "associated_work",
                        "temperature",
                        "pressure",
                    ]
                )["density"]
                .mean()
                .reset_index()
            )

        stackDF = _mask_df(stackDF)

        if sta["combMCMD"]:
            groupREList = ["molecule", "temperature", "forcefield", "pressure"]
        else:
            groupREList = [
                "MCorMD",
                "molecule",
                "temperature",
                "forcefield",
                "pressure",
            ]
        stackDF["Relative_Error"] = stackDF.groupby(groupREList)[
            "density"
        ].transform(_calculate_relative_error)
        # print(stackDF.loc[stackDF["Relative_Error"].abs() > 3])

        _transform_ff_names(stackDF)
        x_label = sta["splitby"]  #'molecule', "forcefield", "associated_work"
        y_label = "Relative_Error"
        color_labels, xtick_labelsList, color_map = _sort_plotting_labels(
            stackDF, x_label, y_label
        )

        if sta["plot_type"] == "scatter":
            _plot_scatter(
                stackDF, fig, ax, color_labels, color_map, x_label, y_label
            )
        elif sta["plot_type"] == "violin":
            _plot_violin(
                stackDF, fig, ax, color_labels, color_map, x_label, y_label
            )
        else:
            return
            # color_labels = list(Counter(df[x_label].to_numpy()).keys())

        _general_plot_formatting(
            stackDF,
            ax,
            sta["splitby"],
            sta["figdir"],
            sta["combMCMD"],
            xtick_labelsList,
            color_labels,
            sta["average_replicas"],
            ymax=ymax,
            rotation=sta.get("rotation", 0),
            wrap_labels=sta.get("wrap_labels", True),
        )  # format all plots
    # remove extra labels
    axs[0].set_ylabel("")
    axs[2].set_ylabel("")

    fig_names = []
    for sta in stack:
        fig_names.append("-".join(list(map(str, sta.values()))))
    fname = f"figures/Stacked_Figures/{'_'.join(fig_names)}.pdf"
    print(f"Plotting {fname}")
    fig.tight_layout()
    fig.savefig(fname, dpi=500, bbox_inches="tight")


def stack_plot_multix(stack):
    """Stack multiple plots on each other with differing xaxis labels on each."""
    # Initialize Figure
    fig, axs = _plot_stack_params_multix(stack)

    # Iterate through stack to generate plot
    for ax, sta in zip(axs, stack):
        ax.axvline(
            x=sta["seperator_line_pos"],
            color="black",
            alpha=0.75,
            linestyle="--",
            zorder=-2,
        )
        # generate filteredDF
        if sta["figdir"] in ["MoSDeF", "RR"]:
            df = sta["df"]
            stackDF = df.loc[df["associated_work"] == sta["figdir"]].copy()
        else:
            df = sta["df"]
            stackDF = df.copy()

        if sta["average_replicas"]:
            stackDF = copy.deepcopy(
                stackDF.groupby(
                    [
                        "molecule",
                        "forcefield",
                        "engine",
                        "MCorMD",
                        "associated_work",
                        "temperature",
                        "pressure",
                    ]
                )["density"]
                .mean()
                .reset_index()
            )

        stackDF = _mask_df(stackDF)

        if sta["combMCMD"]:
            groupREList = ["molecule", "temperature", "forcefield", "pressure"]
        else:
            groupREList = [
                "MCorMD",
                "molecule",
                "temperature",
                "forcefield",
                "pressure",
            ]
        stackDF["Relative_Error"] = stackDF.groupby(groupREList)[
            "density"
        ].transform(_calculate_relative_error)
        print(
            "Number of larger values: ",
            len(stackDF.loc[stackDF["Relative_Error"].abs() > 0.70].index),
        )
        print(
            "Number of smaller values: ",
            len(stackDF.loc[stackDF["Relative_Error"].abs() < 0.70].index),
        )
        # print(stackDF.loc[stackDF["Relative_Error"].abs() > 0.50])

        _transform_ff_names(stackDF)
        x_label = sta["splitby"]  #'molecule', "forcefield", "associated_work"
        y_label = "Relative_Error"
        color_labels, xtick_labelsList, color_map = _sort_plotting_labels(
            stackDF, x_label, y_label
        )

        if sta["plot_type"] == "scatter":
            _plot_scatter(
                stackDF, fig, ax, color_labels, color_map, x_label, y_label
            )
        elif sta["plot_type"] == "violin":
            _plot_violin(
                stackDF, fig, ax, color_labels, color_map, x_label, y_label
            )
        else:
            return
            color_labels = list(Counter(df[x_label].to_numpy()).keys())

        _general_plot_formatting(
            stackDF,
            ax,
            sta["splitby"],
            sta["figdir"],
            sta["combMCMD"],
            xtick_labelsList,
            color_labels,
            sta["average_replicas"],
            ymax=sta["ymax"],
        )  # format all plots
    axs[0].set_ylabel("")
    axs[2].set_ylabel("")

    fname = f"figures/Stacked_Figures/{'molecule_stack'}-{sta['figdir']}.pdf"
    print(f"Plotting {fname}")
    fig.tight_layout()
    fontsize = 20
    axs[2].text(3.65, 5.8, "MoSDeF", fontdict={"size": fontsize})
    axs[2].text(3.65, 29.7, "MoSDeF", fontdict={"size": fontsize})
    axs[2].text(3.65, 53.8, "MoSDeF", fontdict={"size": fontsize})
    axs[2].text(3.1, 5.8, "RR", fontdict={"size": fontsize})
    axs[2].text(3.1, 29.7, "RR", fontdict={"size": fontsize})
    axs[2].text(0.8, 53.8, "RR", fontdict={"size": fontsize})
    fig.savefig(fname, dpi=500, bbox_inches="tight")


def looper_for_plotting_single_data(df):
    """Seperate out plotting data from source specific dataframes."""
    for filter_method in ["MoSDeF", "RR", None]:
        if filter_method:
            df = copy.deepcopy(densityDF)
            df = df.loc[df["associated_work"] == filter_method]
        for use_replica_avg in (True, False):
            for splitby in ["ff-mol", "forcefield", "associated_work"]:
                for combMCMDMCMD in [True, False]:
                    for plot_type in ["scatter", "violin", "box"]:
                        if filter_method is None:
                            single_plot_comparisons(
                                df.copy(),
                                splitby=splitby,
                                combMCMD=combMCMD,
                                plot_type=plot_type,
                                average_replicas=use_replica_avg,
                                figdir="combined",
                            )
                        else:
                            single_plot_comparisons(
                                df.copy(),
                                splitby=splitby,
                                combMCMD=combMCMDMCMD,
                                plot_type=plot_type,
                                average_replicas=use_replica_avg,
                                figdir=filter_method,
                            )
                        df = copy.deepcopy(densityDF)


def plot_rr_mosdef_mosdef_avg(fullDF):
    """Generate the overall averages for intrawork sources and get a total distribution."""
    ylabelfs = 20

    splitby = "associated_work"
    combMCMD = False
    plot_type = "violin"
    average_replicas = [False, False, True]
    figdir = "combined"
    df = fullDF.copy()
    df = _mask_df(df)
    groupREList = [
        "MCorMD",
        "molecule",
        "temperature",
        "forcefield",
        "pressure",
    ]
    df["Relative_Error"] = df.groupby(groupREList)["density"].transform(
        _calculate_relative_error
    )
    _transform_ff_names(df)
    x_label = splitby  #'molecule', "forcefield", "associated_work"
    y_label = "Relative_Error"

    # labels should be Round Robin, MoSDeF, MoSDeF-Averaged
    # color_labels, xtick_labelsList, color_map = _sort_plotting_labels(df, x_label, y_label)
    color_labels = ["RR", "MoSDeF", "MoSDeF-Averaged"]
    dfs = list(
        map(
            lambda x: df.loc[df["associated_work"] == x].copy(),
            color_labels[:2],
        )
    )

    # generate final plotting df
    df3 = fullDF.copy()
    df = (
        df.groupby(
            [
                "molecule",
                "forcefield",
                "engine",
                "MCorMD",
                "associated_work",
                "temperature",
                "pressure",
            ]
        )["density"]
        .mean()
        .reset_index()
    )
    df = _mask_df(df)
    groupREList = [
        "MCorMD",
        "molecule",
        "temperature",
        "forcefield",
        "pressure",
    ]
    df["Relative_Error"] = df.groupby(groupREList)["density"].transform(
        _calculate_relative_error
    )
    dfs.append(df.loc[df["associated_work"] == "MoSDeF"])
    fig, ax = _plot_params(color_labels)

    if plot_type == "scatter":
        _plot_scatter(
            stackDF, fig, ax, color_labels, color_map, x_label, y_label
        )
    elif plot_type == "violin":
        x_ticks = np.arange(len(color_labels))
        for i, iDF in enumerate(dfs):
            series = iDF[y_label]
            ydata = series[series.notnull()]
            if len(x_ticks) > 1:
                spread = (x_ticks[1] - x_ticks[0]) * 0.25
            else:
                spread = 0.25
            print(
                f"Maximum deviation is: {max(abs(ydata)):.2f} for grouping {color_labels[i]}"
            )
            offsets = np.linspace(-1 * spread, spread, len(ydata))
            xdata = np.full(len(ydata), x_ticks[i]) + offsets
            ax.violinplot([ydata], positions=[i], showextrema=True)

    xtick_labelsList = ["Round Robin", "MoSDeF", "MoSDeF-Averaged"]
    _general_plot_formatting(
        dfs[0],
        ax,
        splitby,
        figdir,
        combMCMD,
        xtick_labelsList,
        color_labels,
        average_replicas,
        ymax=8,
    )  # format all plots
    ax.set_ylabel(r"$100\times\Delta\rho / \rho$", size=ylabelfs)

    fname = f"figures/Stacked_Figures/{'Grouped-Comparison'}-{plot_type}.pdf"
    print(f"Plotting {fname}")
    fig.tight_layout()
    fig.savefig(fname, dpi=500, bbox_inches="tight")


def print_errors_for_text(df):
    """Print out relative errors for discussion and conclusion in the manuscript."""
    df = _mask_df(df)
    # std error FF RR
    df_RR = df.loc[(df["associated_work"] == "RR")].copy()
    groupREList = [
        "MCorMD",
        "molecule",
        "temperature",
        "forcefield",
        "pressure",
    ]
    df_RR["Relative_Error"] = df_RR.groupby(groupREList)["density"].transform(
        _calculate_relative_error
    )
    RR_sim_std = df_RR["Relative_Error"].std()
    print(
        f"Round Robin STD Error is {df_RR['Relative_Error'].std()}, Max Error is {df_RR['Relative_Error'].abs().max()}"
    )

    # std error FF MoSDeF
    df_mosdef = df.loc[(df["associated_work"] == "MoSDeF")].copy()
    groupREList = [
        "MCorMD",
        "molecule",
        "temperature",
        "forcefield",
        "pressure",
    ]
    df_mosdef["Relative_Error"] = df_mosdef.groupby(groupREList)[
        "density"
    ].transform(_calculate_relative_error)
    mosdef_sim_std = df_mosdef["Relative_Error"].std()
    print(
        f"MoSDeF STD Error is {df_mosdef['Relative_Error'].std()}, Max Error is {df_mosdef['Relative_Error'].abs().max()}"
    )

    # std error FF MoSDeF average of 16
    groupREList = [
        "MCorMD",
        "molecule",
        "temperature",
        "forcefield",
        "pressure",
    ]
    df_mosdef_avg = df.loc[(df["associated_work"] == "MoSDeF")].copy()
    df_mosdef_avg = copy.deepcopy(
        df_mosdef_avg.groupby(
            [
                "molecule",
                "forcefield",
                "engine",
                "MCorMD",
                "associated_work",
                "temperature",
                "pressure",
            ]
        )["density"]
        .mean()
        .reset_index()
    )
    df_mosdef_avg["Relative_Error"] = df_mosdef_avg.groupby(groupREList)[
        "density"
    ].transform(_calculate_relative_error)
    mosdef_avg_std = df_mosdef_avg["Relative_Error"].std()
    print(
        f"MoSDeF Best Practices STD Error is {df_mosdef_avg['Relative_Error'].std()}, Max Error is {df_mosdef_avg['Relative_Error'].abs().max()}"
    )

    print(f"Comparable MoSDeF Decrease: {RR_sim_std/mosdef_sim_std}")
    print(f"Best Practices MoSDeF Decrease: {RR_sim_std/mosdef_avg_std}")
    mean = df_mosdef_avg["Relative_Error"].mean()
    length_RE = len(df_mosdef_avg["Relative_Error"])
    print(
        f"MoSDeF CI for all errors: {mean} +- {mosdef_avg_std/np.sqrt(length_RE)*1.96}"
    )


if __name__ == "__main__":
    mosdefDF = pd.read_csv("csvs/job_density_data.csv", index_col=0)
    mosdefDF.insert(
        len(mosdefDF.columns),
        "pressure",
        np.full(len(mosdefDF.index), 0.101325),
    )

    hasse_dfList = _load_all_rr_data()

    # 4 combine data into one dataframe
    densityDF = mosdefDF.loc[
        mosdefDF["molecule"] != "pentaneUA-flexible_bonds"
    ].copy()
    densityDF = pd.concat([densityDF, *hasse_dfList], ignore_index=True)
    df = copy.deepcopy(densityDF)

    # generate individual scatter and violin plots
    # looper_for_plotting_single_data(df)

    # generate stacked plots

    stack1 = [
        {
            "figdir": "MoSDeF",
            "splitby": "ff-mol",
            "combMCMD": True,
            "plot_type": "violin",
            "average_replicas": False,
        },
        {
            "figdir": "MoSDeF",
            "splitby": "ff-mol",
            "combMCMD": False,
            "plot_type": "violin",
            "average_replicas": False,
        },
        {
            "figdir": "MoSDeF",
            "splitby": "ff-mol",
            "combMCMD": False,
            "plot_type": "violin",
            "average_replicas": True,
        },
    ]

    stack_plot_comparisons(df, stack1, ymax=0.8)

    stack2 = [
        {
            "figdir": "combined",
            "splitby": "ff-mol",
            "combMCMD": True,
            "plot_type": "scatter",
            "average_replicas": False,
            "rotation": 45,
            "wrap_labels": False,
        },
        {
            "figdir": "combined",
            "splitby": "ff-mol",
            "combMCMD": False,
            "plot_type": "scatter",
            "average_replicas": False,
            "rotation": 45,
            "wrap_labels": False,
        },
    ]

    # stack_plot_comparisons(df, stack2) # try other plotting

    stack3 = [
        {
            "figdir": "combined",
            "splitby": "ff-mol",
            "combMCMD": True,
            "plot_type": "violin",
            "average_replicas": False,
            "rotation": 45,
            "wrap_labels": False,
        },
        {
            "figdir": "combined",
            "splitby": "ff-mol",
            "combMCMD": False,
            "plot_type": "violin",
            "average_replicas": False,
            "rotation": 45,
            "wrap_labels": False,
        },
    ]

    # stack_plot_comparisons(df, stack3)

    # compare TraPPE -> all
    df1 = copy.deepcopy(df).loc[
        df["forcefield"].isin(["TraPPE", "trappe-ua", "benzene-ua"])
    ]
    df2 = copy.deepcopy(df).loc[df["forcefield"].isin(["OPLS", "spce"])]
    df3 = copy.deepcopy(df).loc[df["forcefield"].isin(["OPLSAMBER", "oplsaa"])]
    stack_mol_by_mol_basis = [  # Add text above MoSDeF and RR sections
        {
            "df": df1,
            "figdir": "combined",
            "splitby": "ff-mol",
            "combMCMD": False,
            "plot_type": "violin",
            "average_replicas": False,
            "ymax": 2,
            "seperator_line_pos": 3.5,
        },
        {
            "df": df2,
            "figdir": "combined",
            "splitby": "ff-mol",
            "combMCMD": False,
            "plot_type": "violin",
            "average_replicas": False,
            "ymax": 2,
            "seperator_line_pos": 3.5,
        },
        {
            "df": df3,
            "figdir": "combined",
            "splitby": "ff-mol",
            "combMCMD": False,
            "plot_type": "violin",
            "average_replicas": False,
            "ymax": 8,
            "seperator_line_pos": 3.5,
        },
    ]

    # Compare TraPPE -> TraPPE
    # compare OPLS -> SPC/E
    # compare OPLS AMBER -> OPLS Ethanol
    stack_plot_multix(stack_mol_by_mol_basis)

    plot_rr_mosdef_mosdef_avg(densityDF.copy())

    print_errors_for_text(densityDF.copy())
