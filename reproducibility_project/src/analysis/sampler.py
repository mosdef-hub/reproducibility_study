"""Use the pymbar package to perform decorrelated equilibration sampling."""

import numpy as np
from pymbar import timeseries

from reproducibility_project.src.analysis.equilibration import is_equilibrated


def sample_job(job, variable="potential_energy", threshold_fraction=0.75, threshold_neff=100):
    """Use the timeseries module from pymbar to perform statistical sampling.

    The start, end and decorrleated step size of the production region are
    added to the job document.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The Job object.
    variable : str; default "potential_energy"
        The variable to be used in sampling.
    threshold_fraction : float, optional, default=0.75
        Fraction of data expected to be equilibrated.
    threshold_neff : int, optional, default=100
        Minimum amount of uncorrelated samples to be considered equilibrated
    """
    try:
        job.doc["sampling_results"]
    except KeyError:
        job.doc["sampling_results"] = {}

    data = np.genfromtxt(job.fn("log.txt"), names=True)[variable]
    start, stop, step, Neff = _decorr_sampling(data, threshold_fraction=threshold_fraction, threshold_neff=threshold_neff)
    job.doc["sampling_results"][variable] = (range(start, stop, step), Neff)


def _decorr_sampling(data, threshold_fraction=0.75, threshold_neff=100):
    """Use the timeseries module from pymbar to perform statistical sampling.

    Parameters
    ----------
    data : numpy.Array
        1-D time dependent data to check for equilibration
    threshold_fraction : float, default=0.75
        Fraction of data expected to be equilibrated.
    threshold_neff : int, default=100
        Minimum amount of uncorrelated samples to be considered equilibrated
    """
    is_equil, prod_start, ineff, Neff = is_equilibrated(
        data, threshold_fraction=threshold_fraction, threshold_neff=threshold_neff, nskip=1
    )
    if is_equil:
        uncorr_indices = timeseries.subsampleCorrelatedData(
            data[prod_start:], g=ineff, conservative=True
        )
        return (
            uncorr_indices.start + prod_start,
            uncorr_indices.stop + prod_start,
            uncorr_indices.step,
            Neff,
        )
    else:
        raise ValueError(
            "Property does not have requisite threshold of production data expected."
            "More production data is needed, or the threshold needs to be lowered."
            "See project.src.analysis.equilibration.is_equilibrated for more information."
        )
