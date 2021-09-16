"""Timeseries and pyMBAR related methods."""
import pathlib
from typing import List

import numpy as np
import numpy.typing as npt
import pandas as pd
from pymbar import timeseries
from signac.contrib.job import Job


def is_equilibrated(
    a_t: npt.ArrayLike,
    threshold_fraction: float = 0.8,
    threshold_neff: int = 100,
    nskip: int = 1,
) -> List:
    """Check if a dataset is equilibrated based on a fraction of equil data.

    Using `pymbar.timeseries` module, check if a timeseries dataset has enough
    equilibrated data based on two threshold values. The threshold_fraction
    value translates to the fraction of total data from the dataset 'a_t' that
    can be thought of as being in the 'production' region. The threshold_neff
    is the minimum amount of effectively uncorrelated samples to have in a_t to
    consider it equilibrated.

    The `pymbar.timeseries` module returns the starting index of the
    'production' region from 'a_t'. The fraction of 'production' data is
    then compared to the threshold value. If the fraction of 'production' data
    is >= threshold fraction this will return a list of
    [True, t0, g, Neff] and [False, None, None, None] otherwise.

    Parameters
    ----------
    a_t : numpy.typing.Arraylike
        1-D time dependent data to check for equilibration.
    threshold_fraction : float, optional, default=0.8
        Fraction of data expected to be equilibrated.
    threshold_neff : int, optional, default=100
        Minimum amount of effectively correlated samples to consider a_t
        'equilibrated'.
    nskip : int, optional, default=1
        Since the statistical inefficiency is computed for every time origin
        in a call to timeseries.detectEquilibration, for larger datasets
        (> few hundred), increasing nskip might speed this up, while
        discarding more data.
    """
    if threshold_fraction < 0.0 or threshold_fraction > 1.0:
        raise ValueError(
            f"Passed 'threshold_fraction' value: {threshold_fraction}, "
            "expected value between 0.0-1.0."
        )

    threshold_neff = int(threshold_neff)
    if threshold_neff < 1:
        raise ValueError(
            f"Passed 'threshold_neff' value: {threshold_neff}, expected value "
            "1 or greater."
        )

    [t0, g, Neff] = timeseries.detectEquilibration(a_t, nskip=nskip)
    frac_equilibrated = 1.0 - (t0 / np.shape(a_t)[0])

    if (frac_equilibrated >= threshold_fraction) and (Neff >= threshold_neff):
        return [True, t0, g, Neff]
    else:
        return [False, None, None, None]


def trim_non_equilibrated(
    a_t: npt.ArrayLike,
    threshold_fraction: float = 0.75,
    threshold_neff: int = 100,
    nskip: int = 1,
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
    threshold_fraction : float, optional, default=0.75
        Fraction of data expected to be equilibrated.
    threshold_neff : int, optional, default=100
        Minimum amount of uncorrelated samples.
    nskip : int, optional, default=1
        Since the statistical inefficiency is computed for every time origin
        in a call to timeseries.detectEquilibration, for larger datasets
        (> few hundred), increasing nskip might speed this up, while
        discarding more data.
    """
    [truth, t0, g, Neff] = is_equilibrated(
        a_t,
        threshold_fraction=threshold_fraction,
        threshold_neff=threshold_neff,
        nskip=nskip,
    )
    if not truth:
        raise ValueError(
            f"Data with a threshold_fraction of {threshold_fraction} and "
            f"threshold_neff {threshold_neff} is not equilibrated!"
        )

    return [a_t[t0:], t0, g, Neff]


def plot_job_property_with_t0(
    job: Job,
    filename: str,
    property_name: str,
    log_filename: str = "log.txt",
    title: str = None,
    vline_scale: float = 1.1,
    threshold_fraction: float = 0.0,
    threshold_neff: int = 1,
    overwrite: bool = False,
    data_plt_kwargs: dict = None,
    vline_plt_kwargs: dict = None,
) -> None:
    """Plot data with a vertical line at beginning of equilibration for a specific job and property.

    Parameters
    ----------
    job : signac.contrib.job.Job, required
        The signac job to access the necessary data files.
    filename : str, required
        The name of the output image.
        Only the name of the file and extension is expected, the location will
        be within the job.
    property_name : str, required
        The name of the property to plot.
    log_filename : str, default "log.txt"
        The relative path (from the job directory) to the log file name to read.
    title : str, optional, default = Property
        Title of the plot
    vline_scale : float, optional, default=1.1
        Scale the min and max components of the vertical line.
    threshold_fraction : float, optional, default=0.0
        Fraction of data expected to be equilibrated.
    threshold_neff : int, optional, default=1
        Minimum amount of uncorrelated samples.
    overwrite : bool, optional, default=False
        Do not write to filename if a file already exists with the same name.
        Set to True to overwrite exisiting files.
    data_plt_kwargs : dict, optional, default={}
        Pass in a dictionary of keyword arguments to plot the data.
    vline_plt_kwargs : dict, optional, default={}
        Pass in a dictionary of keyword arguments for the vertical line
        denoting t0.
    """
    from reproducibility_project.src.utils.plotting import (
        plot_data_with_t0_line,
    )

    fname = pathlib.Path(filename)
    fname = fname.name
    a_t = pd.read_csv(
        job.fn(log_filename),
        delim_whitespace=True,
        header=0,
    )
    if data_plt_kwargs is None:
        data_plt_kwargs = dict()
    if vline_plt_kwargs is None:
        vline_plt_kwargs = dict()
    with job:
        plot_data_with_t0_line(
            filename=fname,
            a_t=a_t[property_name].to_numpy(),
            vline_scale=vline_scale,
            title=title,
            overwrite=overwrite,
            threshold_fraction=threshold_fraction,
            threshold_neff=threshold_neff,
            data_plt_kwargs=data_plt_kwargs,
            vline_plt_kwargs=vline_plt_kwargs,
        )
