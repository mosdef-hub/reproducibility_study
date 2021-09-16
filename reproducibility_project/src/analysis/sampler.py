"""Use the pymbar package to perform decorrelated equilibration sampling."""

import numpy as np
import pandas as pd
import signac
from pymbar import timeseries

from reproducibility_project.src.analysis.equilibration import is_equilibrated


def sample_job(
    job,
    filename="log.txt",
    variable="potential_energy",
    threshold_fraction=0.75,
    threshold_neff=100,
):
    """Use the timeseries module from pymbar to perform statistical sampling.

    The start, end and decorrleated step size of the production region are
    added to the job document.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The Job object.
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
    try:
        job.doc["sampling_results"]
    except KeyError:
        job.doc["sampling_results"] = {}

    data = np.genfromtxt(job.fn(filename), names=True)[variable]
    start, stop, step, Neff = _decorr_sampling(
        data,
        threshold_fraction=threshold_fraction,
        threshold_neff=threshold_neff,
    )
    job.doc["sampling_results"][variable] = {
        "start": start,
        "stop": stop,
        "step": step,
        "Neff": Neff,
    }


def write_subsampled_values(
    job: signac.contrib.project.Job,
    property: str,
    property_filename: str = "log.txt",
    overwrite: bool = False,
) -> None:
    """Write out subsampled values to a Job's 'data' document."""
    if not isinstance(job, signac.contrib.project.Job):
        raise TypeError(
            f"Expected input 'job' of type signac.contrib.project.Job, was provided: {type(job)}"
        )

    if property is None or property == "":
        raise ValueError(
            f"Expected 'property' to be a name of a property, was provided {property}."
        )

    if (
        not overwrite
        and job.data.get(f"subsamples/{property}", None) is not None
    ):
        raise ValueError(
            f"Attempting to overwrite already existing data for property: {property}, set overwrite=True to do this."
        )

    sampling_dict = job.doc["sampling_results"][f"{property}"]
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
        property_subsamples = df[f"{property}"].to_numpy()[indices]
        job.data[f"subsamples/{property}"] = property_subsamples


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
