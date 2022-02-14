"""Additional utilities for working with rigid bodies."""

import numpy as np


def moit(points, masses, center=np.zeros(3)):
    """Calculate moment of inertia tensor (moit) for rigid bodies.

    Assumes rigid body center is at origin unless center is provided.
    Only calculates diagonal elements.

    Parameters
    ----------
    points : numpy.ndarray (N,3)
        x, y, and z coordinates of the rigid body constituent particles
    masses : numpy.ndarray (N,)
        Masses of the constituent particles
    center : numpy.ndarray (3,), default np.array([0,0,0])
        x, y, and z coordinates of the rigid body center

    Returns
    -------
    numpy.ndarray (3,)
        moment of inertia tensor for the rigid body center
    """
    points -= center
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    I_xx = np.sum((y**2 + z**2) * masses)
    I_yy = np.sum((x**2 + z**2) * masses)
    I_zz = np.sum((x**2 + y**2) * masses)
    return np.array((I_xx, I_yy, I_zz))
