"""Use the pymbar package to perform decorrelated equilibration sampling."""

from typing import List
from warnings import warn

import numpy as np
import pandas as pd
import signac
import pymbar
from pymbar import timeseries

from reproducibility_project.src.analysis.equilibration import is_equilibrated


def sample_job(
    job: signac.contrib.job.Job,
    ensemble: str,
    filename: str = "log.txt",
    variable: str = "potential_energy",
    threshold_fraction: float = 0.75,
    threshold_neff: int = 100,
    strict: bool = False,
    monte_carlo_override: bool = False,
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
    strict : bool, default=False
        If strict, require both threshold_fraction and threshold_neff to be
        true to evaluate as 'equilibrated'.
    monte_carlo_override : bool, optional, default=False
        Consider the entire data set passed in as production data that is fully equilibrated.
    """
    doc_name = f"{ensemble}/sampling_results"

    data = np.genfromtxt(job.fn(filename), names=True)[variable]
    data_shape = data.shape
    if not monte_carlo_override:
        start, stop, step, Neff = _decorr_sampling(
            data,
            threshold_fraction=threshold_fraction,
            threshold_neff=threshold_neff,
            strict=strict,
        )
    else:
        start = 0
        stop = data_shape[0] - 1
        step = 1
        Neff = data_shape[0]

    try:
        job.doc[doc_name]
    except KeyError:
        job.doc[doc_name] = {}

    if start is not None:
        job.doc[doc_name][variable] = {
            "start": start,
            "stop": stop,
            "step": step,
            "Neff": Neff,
        }
    else:
        warn(f"JOB_ID: {job.id}\nProperty {variable} is not equilibrated.")


def get_subsampled_values(
    job: signac.contrib.project.Job,
    prop: str,
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
    prop : str, required
        The property of interest to write out the subsampled data.
    ensemble : str, required
        The ensemble that the data was sampled from.
    property_filename : str, optional, default="log-npt.txt"
        The filename to sample the data from.

    Examples
    --------
    >>> arr = write_subsampled_values(job, prop="potential_energy",
                                ensemble="npt",
                                property_filename="log-npt.txt")
    >>> assert isinstance(arr, np.ndarray)
    """
    if not isinstance(job, signac.contrib.project.Job):
        raise TypeError(
            f"Expected input 'job' of type signac.contrib.project.Job, was provided: {type(job)}"
        )

    if prop is None or prop == "":
        raise ValueError(
            f"Expected 'prop' to be a name of a property, was provided {prop}."
        )

    sampling_dict = job.doc[f"{ensemble}/sampling_results"][f"{prop}"]
    start = sampling_dict["start"]
    stop = sampling_dict["stop"]
    step = sampling_dict["step"]
    indices = [idx for idx in range(start, stop, step)]

    if not job.isfile(f"{property_filename}"):
        raise FileNotFoundError(
            f"File {property_filename} does not exist in {job}'s workspace."
        )

    with job:
        with open(f"{property_filename}", "r") as f:
            line1 = f.readline()
            df = pd.read_csv(f, delim_whitespace=True, names=line1.replace('#', '').split())
        return df[f"{prop}"].to_numpy()[indices]


def _decorr_sampling(
    data, threshold_fraction=0.75, threshold_neff=100, strict=False
):
    """Use the timeseries module from pymbar to perform statistical sampling.

    Parameters
    ----------
    data : numpy.Array
        1-D time dependent data to check for equilibration
    threshold_fraction : float, default=0.75
        Fraction of data expected to be equilibrated.
    threshold_neff : int, default=100
        Minimum amount of uncorrelated samples to be considered equilibrated
    strict : bool, default=False
        If strict, require both threshold_fraction and threshold_neff to be
        true to evaluate as 'equilibrated'.
    """
    is_equil, prod_start, ineff, Neff = is_equilibrated(
        data,
        threshold_fraction=threshold_fraction,
        threshold_neff=threshold_neff,
        strict=strict,
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
        warn(
            "Property does not have requisite threshold of production data "
            "expected. More production data is needed, or the threshold needs "
            "to be lowered. See project.src.analysis.equilibration.is_equilibrated"
            " for more information."
        )
        return (None, None, None, None)


def get_decorr_samples_using_max_t0(
    job: signac.contrib.Project.Job,
    ensemble: str,
    property_filename: str,
    prop: str,
    threshold_fraction: float = 0.75,
    threshold_neff: int = 100,
    is_monte_carlo: bool = False,
) -> List[float]:
    """Return the subsamples of data according to maximum t0."""
    t0 = job.doc.get(f"{ensemble}/max_t0")

    if t0 is None:
        raise ValueError(
            "Max t0 has not been calculated yet, refer to project-analysis.py"
        )

    with job:
        with open(f"{property_filename}", "r") as f:
            line1 = f.readline()
            df = pd.read_csv(f, delim_whitespace=True, names=line1.replace('#', '').split())
            a_t = df[f"{prop}"].to_numpy()[t0:]
            if is_monte_carlo:
                uncorr_indices = [val for val in range(t0, len(a_t))]
            else:
                try:
                    uncorr_indices = timeseries.subsampleCorrelatedData(
                        A_t=a_t,
                        conservative=True,
                    )
                except pymbar.utils.ParameterError as e:
                    warn(f"This is a captured ParameterError from pymbar {e}, most likely due to zeroes being passed for the values. Returning the unsampled data")
                    uncorr_indices = [i for i in range(t0, len(a_t))]

        return a_t[uncorr_indices]
