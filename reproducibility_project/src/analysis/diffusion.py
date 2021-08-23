"""Facilitate calculation of the diffusion coefficient from simulation data."""

import gsd.hoomd
import matplotlib.pyplot as plt
import numpy as np
from freud.box import Box
from freud.msd import MSD


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


def _gsd_msd(gsdfile, skip, stride=1):
    """Compute the MSD given a GSD file."""
    with gsd.hoomd.open(gsdfile) as trajectory:
        boxes = []
        images = []
        positions = []
        for frame in trajectory[skip::stride]:
            images.append(frame.particles.image)
            boxes.append(frame.configuration.box)
            positions.append(frame.particles.position)

    # msd requires that the box size does not change.
    assert all(all(i == boxes[0]) for i in boxes)
    msd = MSD(Box.from_box(boxes[0]))
    msd.compute(positions, images=images)
    return msd
