"""Use the pymbar package to perform decorrelated equilibration sampling."""

import numpy as np
import pandas as pd
import signac
from pymbar import timeseries

from reproducibility_project.src.analysis.equilibration import is_equilibrated


def sample_job(
    job: signac.contrib.job.Job,
    ensemble: str,
    filename: str = "log.txt",
    variable: str = "potential_energy",
    threshold_fraction: float = 0.75,
    threshold_neff: int = 100,
):
    """Use the timeseries module from pymbar to perform statistical sampling.

    The start, end and decorrleated step size of the production region are
    added to the job document.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The Job object.
    ensemble : str
        The ensemble of interest, affects the name of the sampled values in
        the job.doc
    filename : str, default "log.txt"
        The relative path (from the job directory) to the log file to be
        analyzed.
    variable : str; default "potential_energy"
        The variable to be used in sampling.
    threshold_fraction : float, optional, default=0.75
        Fraction of data expected to be equilibrated.
    threshold_neff : int, optional, default=100
        Minimum amount of uncorrelated samples to be considered equilibrated
    """
    doc_name = f"{ensemble}/sampling_results"
    try:
        job.doc[doc_name]
    except KeyError:
        job.doc[doc_name] = {}

    data = np.genfromtxt(job.fn(filename), names=True)[variable]
    start, stop, step, Neff = _decorr_sampling(
        data,
        threshold_fraction=threshold_fraction,
        threshold_neff=threshold_neff,
    )
    job.doc[doc_name][variable] = {
        "start": start,
        "stop": stop,
        "step": step,
        "Neff": Neff,
    }


def get_subsampled_values(
    job: signac.contrib.project.Job,
    property: str,
    ensemble: str,
    property_filename: str = "log-npt.txt",
) -> None:
    """Return subsampled values based on the sampling results of sample_job.

    Using the results from `sample_job` in the job document, iterate through
    the property file and return a numpy array of the subsampled data.

    This only writes out the subsampled values, the statistical averaging
    will take place in the `analysis-project.py` file.

    Parameters
    ----------
    job : signac.contrib.project.Job, required
        The signac job to operate on.
    property : str, required
        The property of interest to write out the subsampled data.
    ensemble : str, required
        The ensemble that the data was sampled from.
    property_filename : str, optional, default="log-npt.txt"
        The filename to sample the data from.

    Examples
    --------
    >>> arr = write_subsampled_values(job, property="potential_energy",
                                ensemble="npt",
                                property_filename="log-npt.txt")
    >>> assert isinstance(arr, np.ndarray)
    """
    if not isinstance(job, signac.contrib.project.Job):
        raise TypeError(
            f"Expected input 'job' of type signac.contrib.project.Job, was provided: {type(job)}"
        )

    if property is None or property == "":
        raise ValueError(
            f"Expected 'property' to be a name of a property, was provided {property}."
        )

    sampling_dict = job.doc[f"{ensemble}/sampling_results"][f"{property}"]
    start = sampling_dict["start"]
    stop = sampling_dict["stop"]
    step = sampling_dict["step"]
    indices = [idx for idx in range(start, stop, step)]

    if not job.isfile(f"{property_filename}"):
        raise FileNotFoundError(
            f"File {property_filename} does not exist in {job}'s workspace."
        )

    with job:
        df = pd.read_csv(
            f"{property_filename}", delim_whitespace=True, header=0
        )
        return df[f"{property}"].to_numpy()[indices]


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
        data,
        threshold_fraction=threshold_fraction,
        threshold_neff=threshold_neff,
        nskip=1,
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
            "Property does not have requisite threshold of production data "
            "expected. More production data is needed, or the threshold needs "
            "to be lowered. See project.src.analysis.equilibration.is_equilibrated"
            " for more information."
        )
