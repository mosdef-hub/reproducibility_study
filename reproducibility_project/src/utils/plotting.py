"""Utility methods for plotting data."""
import pathlib
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import scipy

from reproducibility_project.src.analysis.equilibration import is_equilibrated


def plot_data_with_t0_line(
    filename: str,
    a_t: npt.ArrayLike,
    vline_scale: float = 1.1,
    threshold_fraction: float = 0.0,
    threshold_neff: int = 1,
    overwrite: bool = False,
    title: str = None,
    data_plt_kwargs: dict = None,
    vline_plt_kwargs: dict = None,
) -> None:
    """Plot data with a vertical line at beginning of equilibration.

    Parameters
    ----------
    a_t : numpy.typing.ArrayLike
        1-D time dependent data
    vline_scale : float, optional, default=1.1
        Scale the min and max components of the vertical line.
    threshold_fraction : float, optional, default=0.0
        Threshold_fraction to error out on if threshold fraction of data is not equilibrated.
    threshold_neff : int, optional, default=1
        Threshold amount of uncorrelated samples to error out on if threshold amount of data is not equilibrated.
    overwrite : bool, optional, default=False
        Do not write to filename if a file already exists with the same name.
        Set to True to overwrite exisiting files.
    title : str, optional
        Set the title of the plot
    data_plt_kwargs : dict, optional, default={}
        Pass in a dictionary of keyword arguments to plot the data.
    vline_plt_kwargs : dict, optional, default={}
        Pass in a dictionary of keyword arguments for the
    """
    path = pathlib.Path(filename)
    if path.is_file() and not overwrite:
        raise FileExistsError(
            f"Cannot write {path.name}, file already exists at: {path.absolute()}. Set overwrite=True to replace file."
        )

    _, t0, g, Neff = is_equilibrated(
        a_t,
        threshold_fraction=threshold_fraction,
        threshold_neff=threshold_neff,
        nskip=1,
    )

    ymin = np.min(a_t) - np.min(a_t) * (np.abs(1 - vline_scale))
    ymax = np.max(a_t) + np.max(a_t) * (np.abs(1 - vline_scale))

    if data_plt_kwargs is None:
        data_plt_kwargs = {}
    if vline_plt_kwargs is None:
        vline_plt_kwargs = {}
    for key, val in {
        "color": "b",
        "linestyle": "-",
        "label": "Property",
    }.items():
        if data_plt_kwargs.get(key) is None:
            data_plt_kwargs[key] = val

    for key, val in {
        "colors": "r",
        "linestyles": "--",
        "label": "t_0={0:d}\ng={1:.2f}\nNeff={2:d}",
    }.items():
        if vline_plt_kwargs.get(key) is None:
            if key == "label":
                vline_plt_kwargs[key] = val.format(t0, g, int(Neff))
            else:
                vline_plt_kwargs[key] = val

    fig, ax = plt.subplots()

    line1 = ax.plot(a_t, **data_plt_kwargs)
    vline1 = ax.vlines(x=t0, ymin=ymin, ymax=ymax, **vline_plt_kwargs)

    ax.legend(loc="best")

    if title is None:
        plt.suptitle("Property")
    else:
        plt.suptitle(title)

    plt.tight_layout()
    fig.savefig(
        str(path.absolute()),
    )
    plt.close()


def mean_confidence_interval(m, se, confidence=0.95):
    """Calculate mean confidence_interval."""
    se = np.array(se)
    n = 16
    h = se * scipy.stats.t.ppf((1 + confidence) / 2.0, n - 1)
    return h


def wrap_labels(ax, width, break_long_words=False):
    """Wrap labels if string is too long."""
    labels = []
    for label in ax.get_xticklabels():
        text = label.get_text()
        labels.append(
            textwrap.fill(text, width=width, break_long_words=break_long_words)
        )
    ax.set_xticklabels(labels, rotation=0)


# Colors
symbols = {}
symbols["cassandra"] = "o"
symbols["mcccs"] = "^"
symbols["gomc"] = "s"
symbols["gromacs"] = "x"
symbols["hoomd"] = "v"
symbols["lammps-VU"] = "D"
# symbols[ "LAMMPS-UD"] = ">"

colors = {}
colors["cassandra"] = "#009392"  # (23/256, 109/256, 156/256)
colors["mcccs"] = "#39B1B5"  # (194/256, 135/256, 32/256)
colors["gomc"] = "#9CCB86"  # (21/256, 138/256, 106/256)
colors["gromacs"] = "#E9E29C"  # (186/256, 97/256, 26/256)
colors["hoomd"] = "#EEB479"  # (193/256, 130/256, 181/256)
colors["lammps-VU"] = "#E88471"  # (188/256, 146/256, 110/256)
colors["fix"] = "#FF7F00"
colors["flex"] = "#007FFF"
# colors_dict = {all_engine_molecule: all_engine_colors,
#                pentane_fixed: pentane_fixed_colors,
#                pentane_flexible: pentane_flexible_colors}

pretty_names = {
    "mcccsflex": "MCCCS-MN-FlexOH",
    "mcccsfix": "MCCCS-MN-FixOH",
    "lammps-VUflex": "LAMMPS-FlexOH",
    "lammps-VUfix": "LAMMPS-FixOH",
    "cassandra": "Cassandra",
    "mcccs": "MCCCS-MN",
    "gomc": "GOMC",
    "gromacs": "GROMACS",
    "hoomd": "HOOMD-Blue",
    "lammps-VU": "LAMMPS",
}


figsize = (9, 6)
