import os

import foyer
import gsd
import gsd.hoomd
import numpy as np
import pytest
import signac

from project.src import xmls


class BaseTest:
    @pytest.fixture
    def spceff(self, name="spce.xml"):
        abs_path = os.path.dirname(os.path.abspath(xmls.__file__))
        return foyer.Forcefield(forcefield_files=str(abs_path) + "/" + name)

    @pytest.fixture
    def job_gsdfile(self, tmp_job):
        filename = tmp_job.fn("trajectory.gsd")
        create_gsd(filename)
        return tmp_job

    @pytest.fixture
    def tmp_project(self):
        with signac.TemporaryProject(name="test") as p:
            return p

    @pytest.fixture
    def tmp_job(self, tmp_project, statepoint={"a": 0}):
        tmp_project.open_job(statepoint).init()
        for job in tmp_project:
            return job

    def trappe_ua(self):
        return foyer.forcefields.load_TRAPPE_UA()


def create_frame_random(i, bonds, n=50, L=10, seed=42):
    np.random.seed(seed)
    s = gsd.hoomd.Snapshot()
    s.configuration.step = i
    s.particles.N = n
    # First half of particles are "A" and the rest are "B"
    half = n - n // 2
    s.particles.types = ["A", "B"]
    s.particles.typeid = np.zeros(n)
    s.particles.typeid[half:] = 1
    s.particles.position = np.random.random(size=(n, 3)) * L - L / 2
    s.configuration.box = [L, L, L, 0, 0, 0]
    if bonds:
        s.bonds.N = 2
        s.bonds.types = ["AB"]

        # Create bonds between A-B particles within d_cut
        box = freud.box.Box.cube(L)
        aq = freud.locality.AABBQuery(box, s.particles.position[:half])
        d_cut = 2
        group = []
        selection = {"num_neighbors": 1}
        for i, j, d in aq.query(s.particles.position[half:], selection):
            if d < d_cut:
                group.append((j, i + half))
        s.bonds.typeid = np.zeros(len(group))
        s.bonds.group = group
    s.particles.image = np.zeros(shape=(n, 3))
    s.validate()
    return s


def create_gsd(filename, bonds=False, system="random"):
    if system == "random":
        with gsd.hoomd.open(name=filename, mode="wb") as f:
            f.extend(
                [create_frame_random(i, bonds=bonds, seed=i) for i in range(10)]
            )
    else:
        with gsd.hoomd.open(name=filename, mode="wb") as f:
            f.extend([create_frame_xstal(i, seed=i) for i in range(10)])
