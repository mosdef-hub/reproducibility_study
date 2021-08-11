"""Facilitate the calculation of the RDF from simulation data."""

import freud
import gsd
import gsd.hoomd
import matplotlib.pyplot as plt
import numpy as np


def gsd_rdf(job, frames=10, stride=1, bins=50, r_min=0.5, r_max=None):
    """Compute the RDF given a GSD file.

    A wrapper for freud's RDF module
    https://freud.readthedocs.io/en/latest/modules/density.html#freud.density.RDF

    Parameters
    ----------
    job : signac.contrib.job.Job
        The Job object.
    frames : int, default 10
        The number of frames from the trajectory to average. Up to and always
        uses the last frame.
    stride : int, default 1
        The step size between frames
    bins : int, default 50
        The number of bins in the RDF histogram.
    r_min : float, default 0.5
        The minimum distance to calculate the RDF.
    r_max : float, default None
        The maximum distance to calculate the RDF. If None is provided, then
        the minimum box length times a factor of 0.45 will be used.

    TODO (units)
    """
    gsdfile = job.fn("trajectory.gsd")

    if r_max is None:
        with gsd.hoomd.open(gsdfile) as trajectory:
            box = trajectory[-1].configuration.box[:3]
        r_max = min(box) * 0.45

    rdf = freud.density.RDF(bins=bins, r_min=r_min, r_max=r_max)

    with gsd.hoomd.open(gsdfile) as trajectory:
        start = -(frames * stride) + stride - 1
        for frame in trajectory[start::stride]:
            rdf.compute(frame, reset=False)

    fig, ax = plt.subplots()
    ax.plot(rdf.bin_centers, rdf.rdf)
    ax.set_xlabel("$r$")
    ax.set_ylabel("$g(r)$")
    ax.set_title("RDF")

    fig.savefig(job.fn("rdf.png"))

    rdf_array = np.vstack((rdf.bin_centers, rdf.rdf)).T
    np.savetxt(job.fn("rdf.txt"), rdf_array)
