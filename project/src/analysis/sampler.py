"""Use the pymbar package to perform decorrelated equilibration sampling."""

import numpy as np
from pymbar import timeseries


def decorr_sampling(job, data_type="potential_energy"):
    """Uses the timeseries module from pymbar to perform statistical sampling.

    The start, end and decorrleated step size of the production region are
    added to the job document.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The Job object.
    data_type : str; default "potential_energy"
        The type of simulation data to be used in sampling
    """
    data = np.genfromtxt(job.fn("log.txt", names=True))[data_type]
    prod_start, ineff, prod_size = timeseries.detectEquilibration(data)
    uncorr_indices = timeseries.subsampleCorrelatedData(
        data[prod_start:], g=ineff, conservative=True
    )
    job.doc["sample_start"] = uncorr_indices.start + prod_start
    job.doc["sample_end"] = uncorr_indices.stop + prod_start
    job.doc["sample_stride"] = uncorr_indices.step
