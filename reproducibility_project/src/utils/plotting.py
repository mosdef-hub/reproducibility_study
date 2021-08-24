"""Utility methods for plotting data."""
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

from reproducibility_project.src.analysis.equilibration import is_equilibrated


def plot_data_with_t0_line(
    filename: str,
    a_t: npt.ArrayLike,
    vline_scale: float = 1.5,
    threshold: float = 0.0,
    overwrite: bool = False,
    title: str = None,
    data_plt_kwargs: dict = {},
    vline_plt_kwargs: dict = {},
) -> None:
    """Plot data with a vertical line at beginning of equilibration.

    Parameters
    ----------
    a_t : numpy.typing.ArrayLike
        1-D time dependent data
    vline_scale : float, optional, default=1.5
        Scale the min and max components of the vertical line.
    threshold : float, optional, default=0.0
        Threshold to error out on if threshold fraction of data is not equilibrated.
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

    _, t0, g, Neff = is_equilibrated(a_t, threshold=threshold, nskip=1)

    ymin = np.min(a_t) * vline_scale
    ymax = np.max(a_t) * vline_scale

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
        "label": f"t_0={t0}\ng={g:.2f}\nNeff={Neff:.2f}",
    }.items():
        if vline_plt_kwargs.get(key) is None:
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
