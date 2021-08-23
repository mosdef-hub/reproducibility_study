"""Timeseries and pyMBAR related methods."""
import pathlib
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from pymbar import timeseries
from unyt import matplotlib_support


def is_equilibrated(
    a_t: npt.ArrayLike, threshold: float = 0.8, nskip: int = 1
) -> List:
    """Check if a dataset is equilibrated based on a fraction of equil data.

    Using `pymbar.timeseries` module, check if a timeseries dataset has enough
    equilibrated data based on a threshold value. The threshold value
    translates to the fraction of total data from the dataset 'a_t' that
    can be thought of as being in the 'production' region.

    The `pymbar.timeseries` module returns the starting index of the
    'production' region from 'a_t'. The fraction of 'production' data is
    then compared to the threshold value. If the fraction of 'production' data
    is >= threshold fraction this will return a list of
    [True, t0, g, Neff] and [False, None, None, None] otherwise.

    Parameters
    ----------
    a_t : numpy.typing.Arraylike
        1-D time dependent data to check for equilibration.
    threshold : float, optional, default=0.8
        Fraction of data expected to be equilibrated.
    nskip : int, optional, default=1
        Since the statistical inefficiency is computed for every time origin
        in a call to timeseries.detectEquilibration, for larger datasets
        (> few hundred), increasing nskip might speed this up, while
        discarding more data.
    """
    if threshold < 0.0 or threshold > 1.0:
        raise ValueError(
            f"Passed 'threshold' value: {threshold}, expected value between 0.0-1.0."
        )

    [t0, g, Neff] = timeseries.detectEquilibration(a_t, nskip=nskip)
    frac_equilibrated = 1.0 - (t0 / np.shape(a_t)[0])

    if frac_equilibrated >= threshold:
        return [True, t0, g, Neff]
    else:
        return [False, None, None, None]


def trim_non_equilibrated(
    a_t: npt.ArrayLike, threshold: float = 0.75, nskip: int = 1
) -> List:
    """Prune timeseries array to just the production data.

    Refer to equilibration.is_equilibrated for addtional information.

    This method returns a list of length 3, where list[0] is the trimmed array,
    list[1] is the index of the original dataset where equilibration begins,
    list[2] is the calculated statistical inefficiency, which can be used
    when subsampling the data using `pymbar.timseries.subsampleCorrelatedData`,
    list[3] is the number of effective uncorrelated data points.

    Refer to https://pymbar.readthedocs.io/en/master/timeseries.html for
    additional information.

    Parameters
    ----------
    a_t : numpy.typing.Arraylike
        1-D time dependent data to check for equilibration.
    threshold : float, optional, default=0.8
        Fraction of data expected to be equilibrated.
    nskip : int, optional, default=1
        Since the statistical inefficiency is computed for every time origin
        in a call to timeseries.detectEquilibration, for larger datasets
        (> few hundred), increasing nskip might speed this up, while
        discarding more data.

    """
    [truth, t0, g, Neff] = is_equilibrated(
        a_t, threshold=threshold, nskip=nskip
    )
    if not truth:
        raise ValueError(
            f"Data with a threshold of {threshold} is not equilibrated!"
        )

    return [a_t[t0:], t0, g, Neff]


def plot_data_with_t0_line(
    filename: str,
    a_t: npt.ArrayLike,
    threshold: float = 0.0,
    overwrite: bool = False,
    plt_kwargs: dict = None,
) -> None:
    """Plot data with a vertical line at beginning of equilibration.

    Parameters
    ----------
    a_t : numpy.typing.ArrayLike
        1-D time dependent data
    threshold : float, optional, default=0.0
        Threshold to error out on if threshold fraction of data is not equilibrated.
    overwrite : bool, optional, default=False
        Do not write to filename if a file already exists with the same name.
        Set to True to overwrite exisiting files.
    plt_kwargs : dict, optional, default=None
        keyword dictionary for matplotlib.pyplot.plot command.
    """
    scale = 1.5
    path = pathlib.Path(filename)
    if path.is_file() and not overwrite:
        raise FileExistsError(
            f"Cannot write {path.name}, file already exists at: {path.absolute()}. Set overwrite=True to replace file."
        )

    _, t0, g, Neff = is_equilibrated(a_t, threshold=threshold, nskip=1)

    ymin = np.min(a_t) * scale
    ymax = np.max(a_t) * scale

    (line1,) = plt.plot(
        a_t,
        "b-",
        label=f"Property",
    )
    vline1 = plt.vlines(
        x=t0,
        ymin=ymin,
        ymax=ymax,
        colors="r",
        linestyles="--",
        label=f"t_0={t0}\ng={g:.2f}\nNeff={Neff:.2f}",
    )
    plt.legend(loc="best")
    plt.savefig(
        str(path.absolute()),
    )
