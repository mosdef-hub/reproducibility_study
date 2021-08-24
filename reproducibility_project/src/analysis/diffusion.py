"""Facilitate calculation of the diffusion coefficient from simulation data."""

import gsd.hoomd
import matplotlib.pyplot as plt
import numpy as np
import unyt as u
from freud.box import Box
from freud.msd import MSD
from scipy import stats


def gsd_msd(job, skip=2, stride=1):
    """Compute the MSD and diffusion coefficient given a Signac Job object.

    The job folder is expected to contain the file "trajectory.gsd" with lengths
    in nanometers.
    This function is a convenience wrapper for freud's MSD module
    https://freud.readthedocs.io/en/latest/modules/msd.html
    The diffusion coefficient is calculated as the slope/6 of a linear fit to
    the MSD following https://doi.org/10.1002/jcc.21939
    After execution, the files "msd.png" and "msd.txt" are created in the job
    folder and the "diffusion_coefficent" (in nm^2/ps) is set in the job.doc.

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
    msd, timesteps = _gsd_msd(job.fn("trajectory.gsd"), skip, stride)

    m, b, r, p, std_err = stats.linregress(timesteps, msd.msd)

    fig, ax = plt.subplots()
    ax.plot(
        timesteps,
        m * timesteps + b,
        label=f"linear fit\ny = {m:.1e}x + {b:.1e}\n(r = {r:.3f})",
    )
    ax.plot(timesteps, msd.msd, label="MSD")
    ax.set_xlabel("$Time (ps)$")
    ax.set_ylabel(r"$MSD (nm^{2})$")
    ax.set_title("MSD")

    fig.savefig(job.fn("msd.png"))

    np.savetxt(job.fn("msd.txt"), msd.msd)

    # calculated according to
    # https://doi.org/10.1002/jcc.21939
    # units nm^2/ps
    job.doc["diffusion_coefficient"] = m / 6
    return msd


def _gsd_msd(gsdfile, skip, stride=1):
    """Compute the MSD given a GSD file."""
    with gsd.hoomd.open(gsdfile) as trajectory:
        boxes = []
        images = []
        positions = []
        timesteps = []
        for frame in trajectory[skip::stride]:
            boxes.append(frame.configuration.box)
            images.append(frame.particles.image)
            positions.append(frame.particles.position)
            timesteps.append(frame.configuration.step)

    if not all(all(i == boxes[0]) for i in boxes):
        raise ValueError(
            "MSD calculation requires that the box size does not change."
        )
    msd = MSD(Box.from_box(boxes[0]))
    msd.compute(positions, images=images)
    return msd, np.array(timesteps)
