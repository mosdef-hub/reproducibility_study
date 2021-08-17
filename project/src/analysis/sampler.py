"""Use the pymbar package to perform decorrelated equilibration sampling."""

import numpy as np
from pymbar import timeseries

from project.src.analysis.equlibration import is_equilibrated


def decorr_sampling(job, data_type="potential_energy"):
    """Use the timeseries module from pymbar to perform statistical sampling.

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
    is_equil, prod_start, ineff = is_equilibrated(data, threshold=0.8, nskip=1)
    if is_equil:
        uncorr_indices = timeseries.subsampleCorrelatedData(
            data[prod_start:], g=ineff, conservative=True
        )
        job.doc["sample_start"] = uncorr_indices.start + prod_start
        job.doc["sample_end"] = uncorr_indices.stop + prod_start
        job.doc["sample_stride"] = uncorr_indices.step
    else:
        raise ValueError(
            "Property does not have requisite threshold of production data expected."
            "More production data is needed, or the threshold needs to be lowered."
            "See project.src.analysis.equilibration.is_equilibrated for more information."
        )
