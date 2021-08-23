"""Facilitate calculation of the diffusion coefficient from simulation data."""

import freud
import gsd.hoomd
import matplotlib.pyplot as plt
import numpy as np


def gsd_msd(job, skip=2, stride=1):
    """Compute the MSD given a Signac Job object.

    The job folder is expected to contain the file "trajectory.gsd" with lengths
    in nanometers.
    This function is a convenience wrapper for freud's MSD module
    https://freud.readthedocs.io/en/latest/modules/msd.html
    After execution, the files "msd.png" and "msd.txt" are created in the job
    folder.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The Job object.
    skip : int, default 2
        The number of frames from the trajectory to skip.
    stride : int, default 1
        The step size between frames.

    Returns
    -------
    freud.msd.MSD
        Computed MSD object
    """
    msd = _gsd_msd(job.fn("trajectory.gsd"), skip, stride)

    fig, ax = plt.subplots()
    ax.plot(msd.msd)
    ax.set_xlabel("$Frame$")
    ax.set_ylabel(r"$MSD (nm^{2})$")  # TODO check unit
    ax.set_title("MSD")

    fig.savefig(job.fn("msd.png"))

    np.savetxt(job.fn("msd.txt"), msd.msd)
    return msd


def _gsd_rdf(gsdfile, frames=10, stride=1, bins=50, r_min=0.5, r_max=None):
    """Compute the RDF given a GSD file."""
    if r_max is None:
        with gsd.hoomd.open(gsdfile) as trajectory:
            box = trajectory[-1].configuration.box[:3]
        r_max = min(box) * 0.45

    rdf = freud.density.RDF(bins=bins, r_min=r_min, r_max=r_max)

    with gsd.hoomd.open(gsdfile) as trajectory:
        start = -(frames * stride) + stride - 1
        for frame in trajectory[start::stride]:
            rdf.compute(frame, reset=False)

    return rdf
