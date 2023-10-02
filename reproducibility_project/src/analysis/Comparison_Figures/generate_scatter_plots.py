"""Analyze csvs in ./csvs in order to create scatter plots of the density differences found in the SI of 10.1021/acs.jctc.7b00489."""

import camelot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

np.random.seed(0)
import matplotlib.transforms as transforms


def _create_scatter_plot(summary_df, molecule, title, savefig=False):
    """Generate summary files for the aggregate project."""
    plt.rcParams.update(
        {
            "axes.titlesize": 20,
            "axes.labelsize": 20,
            "axes.titleweight": "bold",
            "axes.labelweight": "bold",
            "font.weight": "bold",
            "xtick.labelsize": 20,
            "ytick.labelsize": 20,
        }
    )
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(10, 10))
    offset = lambda p: transforms.ScaledTranslation(
        p / 72, 0, plt.gcf().dpi_scale_trans
    )
    trans = plt.gca().transData

    x_label = "T / K"
    y_labels = ["ρ / kg m−3", "Relative_Error"]
    y_label = y_labels[1]

    color_labels = summary_df["Program(Group)"].unique()
    palette = sns.color_palette("colorblind", n_colors=len(color_labels))
    color_map = dict(zip(color_labels, palette))
    offsets = np.linspace(-12, 12, len(color_labels))
    x_ticks = np.arange(len(color_labels))
    for offseti, my_engine in zip(offsets, color_labels):
        series = summary_df.loc[summary_df["Program(Group)"] == my_engine]
        data = series[series.notnull()]
        ax = plt.scatter(
            x=data[x_label] + offseti,
            y=data[y_label],
            color=color_map[my_engine],
            label=my_engine,
            s=np.full(len(data.index), 40),
            edgecolors="black",
        )
    plt.legend(color_labels)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    val_min = summary_df[y_label].min()
    val_max = summary_df[y_label].max()

    # plt.ylim(val_min - .2 * val_min, val_max + .2 * val_max)
    plt.title(title, y=1.0, pad=15)
    plt.tight_layout()
    if savefig:
        plt.savefig(f"figures/{molecule}_summary_scatter.pdf")
    return fig


def generate_REPlot_by_molecule(molecule, savefig=True):
    """Compare errors in density for engines in csvs."""
    if "AMBER" in molecule:
        return  # only look at molecules with opls and trappe forcefield
    newDF = pd.read_csv(f"csvs/{molecule}_data.csv", index_col=0)
    title = f"SI Data from Hans Hasse for {molecule}"
    x = _create_scatter_plot(newDF, molecule, title, savefig=savefig)


def _generate_csv_from_PDF(molecule):
    """Methods to generate dataframe from PDF supplied for paper SI."""
    page = molecule_pageDictionary[molecule]

    # Read remote pdf into a list of DataFrame
    link = "Hasse_SI.pdf"  # must have access locally for this python file to load it.
    tables = camelot.read_pdf(link, flavor="stream", pages=page)

    if "AMBER" in molecule:
        return

    table = tables[0]
    molecule_nameStr = table.df[1][0]
    if molecule in [
        "OPLS Ethane at 5 MPa",
        "OPLS Propane at 5 MPa",
        "OPLS n-Butane at 5 MPa",
    ]:
        df = table.df[table.df.columns[1:]][2:]
        df.columns = np.arange(len(df.columns))
    elif "AMBER" in molecule:
        df = table.df[table.df.columns[:]][2:]
    else:
        df = table.df[table.df.columns[:]][1:]
    row_start = 1 + np.where(df[0] == "Program(Group)")[0][0]
    headers = df.iloc[row_start - 1]
    newDF = pd.DataFrame(df.values[row_start:], columns=headers)
    for i in range(len(newDF.index)):
        headerInt = int(4 * np.floor(i / 4))
        header = newDF.iloc[headerInt]["Program(Group)"]
        newDF.iloc[i]["Program(Group)"] = header

    # we will change the data type
    # of id column to str by giving
    # the dict to the astype method
    columns = list(newDF.columns[1:])

    for c in columns:
        newDF[c] = pd.to_numeric(newDF[c], errors="coerce")
    newDF["statepoint"] = "temp: " + newDF["T / K"].astype(str)
    newDF["Relative_Error"] = newDF.groupby("statepoint")[
        "ρ / kg m−3"
    ].transform(_calculate_relative_error)
    newDF.to_csv(f"csvs/{molecule}_data.csv")


def _calculate_relative_error(ser):
    """Calculate relative density error for all engines at that statepoint."""
    dataArray = ser.to_numpy()
    mean = np.nanmean(dataArray)
    return (dataArray - mean) / mean * 100


if __name__ == "__main__":
    molecule_pageDictionary = {
        "OPLS Ethane at 5 MPa": "8",
        "OPLS Ethane at 41 MPa": "9",
        "OPLS Ethane at 70 MPa": "10",
        "TraPPE Ethane at 5 MPa": "11",
        "TraPPE Ethane at 41 MPa": "12",
        "TraPPE Ethane at 70 MPa": "13",
        "OPLSAMBER Ethane at 5 MPa": "14,15",
        "OPLSAMBER Ethane at 41 MPa": "15,16",
        "OPLSAMBER Ethane at 70 MPa": "16,17",
        "OPLS Propane at 5 MPa": "18",
        "OPLS Propane at 41 MPa": "19",
        "OPLS Propane at 70 MPa": "20",
        "TraPPE Propane at 5 MPa": "21",
        "TraPPE Propane at 41 MPa": "22",
        "TraPPE Propane at 70 MPa": "23",
        "OPLSAMBER Propane at 5 MPa": "24-25",
        "OPLSAMBER Propane at 41 MPa": "25-26",
        "OPLSAMBER Propane at 70 MPa": "26-27",
        "OPLS n-Butane at 5 MPa": "28",
        "OPLS n-Butane at 41 MPa": "29",
        "OPLS n-Butane at 70 MPa": "30",
        "TraPPE n-Butane at 5 MPa": "31",
        "TraPPE n-Butane at 41 MPa": "32",
        "TraPPE n-Butane at 70 MPa": "33",
        "OPLSAMBER n-Butane at 5 MPa": "34,35",
        "OPLSAMBER n-Butane at 41 MPa": "35,36",
        "OPLSAMBER n-Butane at 70 MPa": "36,37",
    }
    for molecule in molecule_pageDictionary:
        print(f"Generating Plot for: {molecule}")
        print("################################")
        print("################################")
        generate_REPlot_by_molecule(molecule, savefig=True)
