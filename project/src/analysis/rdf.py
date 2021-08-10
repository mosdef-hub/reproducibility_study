"""Facilitate the calculation of the RDF from simulation data."""

import freud
import gsd
import gsd.hoomd
import matplotlib.pyplot as plt
import numpy as np


def gsd_rdf(job):
    """Compute the RDF given a GSD file."""
    gsdfile = job.fn("trajectory.gsd")

    # TODO check that these are reasonable values
    frames = 10
    bins = 50
    r_min = 0.5
    r_max = 2.5

    rdf = freud.density.RDF(bins=bins, r_min=r_min, r_max=r_max)

    with gsd.hoomd.open(gsdfile) as trajectory:
        for frame in trajectory[-frames:-1]:
            rdf.compute(frame, reset=False)

    fig, ax = plt.subplots()
    ax.plot(rdf.bin_centers, rdf.rdf)
    ax.set_xlabel("$r$")
    ax.set_ylabel("$g(r)$")
    ax.set_title("RDF")

    fig.savefig(job.fn("rdf.png"))

    rdf_array = np.vstack((rdf.bin_centers, rdf.rdf)).T
    np.savetxt(job.fn("rdf.txt"), rdf_array)
