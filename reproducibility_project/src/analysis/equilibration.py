"""Timeseries and pyMBAR related methods."""
import pathlib
from typing import List

import numpy as np
import numpy.typing as npt
import pandas as pd
from pymbar import timeseries
from signac.contrib.job import Job


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


def plot_job_property_with_t0(
    job: Job,
    filename: str,
    property_name: str,
    vline_scale: float = 1.5,
    threshold: float = 0.0,
    overwrite: bool = False,
    data_plt_kwargs: dict = {},
    vline_plt_kwargs: dict = {},
) -> None:
    """Plot data with a vertical line at beginning of equilibration for a specifc job and property.

    Parameters
    ----------
    job : signac.contrib.job.Job, required
        The signac job to access the necessary data files.
    filename : str, required
        The name of the output image.
        Only the name of the file and extension is expected, the location will be within the job.
    property : str, required
        The name of the property to plot.
    vline_scale : float, optional, default=1.5
        Scale the min and max components of the vertical line.
    threshold : float, optional, default=0.0
        Threshold to error out on if threshold fraction of data is not equilibrated.
    overwrite : bool, optional, default=False
        Do not write to filename if a file already exists with the same name.
        Set to True to overwrite exisiting files.
    data_plt_kwargs : dict, optional, default={}
        Pass in a dictionary of keyword arguments to plot the data.
    vline_plt_kwargs : dict, optional, default={}
        Pass in a dictionary of keyword arguments for the vertical line denoting t0.
    """
    from reproducibility_project.src.utils.plotting import (
        plot_data_with_t0_line,
    )

    fname = pathlib.Path(filename)
    fname = fname.name
    a_t = pd.read_csv(
        job.fn("log.txt"),
        delim_whitespace=True,
        header=0,
    )
    with job:
        plot_data_with_t0_line(
            filename=fname,
            a_t=a_t[property_name].to_numpy(),
            vline_scale=vline_scale,
            overwrite=overwrite,
            threshold=threshold,
            data_plt_kwargs=data_plt_kwargs,
            vline_plt_kwargs=vline_plt_kwargs,
        )
