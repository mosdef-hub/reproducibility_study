"""Timeseries and pyMBAR related methods."""
from typing import List

import numpy as np
import numpy.typing as npt
from pymbar import timeseries


def is_equilibrated(
    a_t: npt.ArrayLike, threshold: float = 0.8, nskip: int = 1
) -> bool:
    """Check if a dataset is equilibrated based on a fraction of equil data.

    Using `pymbar.timeseries` module, check if a timeseries dataset has enough
    equilibrated data based on a threshold value. The threshold value
    translates to the fraction of total data from the dataset 'a_t' that
    can be thought of as being in the 'production' region.

    The `pymbar.timeseries` module returns the starting index of the
    'production' region from 'a_t'. The fraction of 'production' data is
    then compared to the threshold value. If the fraction of 'production' data
    is >= threshold fraction this will return True, False otherwise.

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

    [t0, _, _] = timeseries.detectEquilibration(a_t, nskip=nskip)
    frac_equilibrated = 1.0 - (t0 / np.shape(a_t)[0])

    if frac_equilibrated >= threshold:
        return True
    else:
        return False


def trim_non_equilibrated(
    a_t: npt.ArrayLike, threshold: float = 0.75, nskip: int = 1
) -> List:
    """Prune timeseries array to just the production data.

    Refer to equilibration.is_equilibrated for addtional information.#!/usr/bin/env python

    This method returns a list of length 3, where list[0] is the trimmed array,
    list[1] is the calculated statistical inefficiency and list[2] is the
    index of the original dataset where equilibration begins, which can be used
    when subsampling the data using `pymbar.timseries.subsampleCorrelatedData`.

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
    if not is_equilibrated(a_t, threshold=threshold, nskip=nskip):
        raise ValueError(
            f"Data with a threshold of {threshold} is not equilibrated!"
        )

    [t0, g, _] = timeseries.detectEquilibration(a_t, nskip=nskip)

    return [a_t[t0:], g, t0]
