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


def create_frame(i, bonds, images, seed=42):
    np.random.seed(seed)
    s = gsd.hoomd.Snapshot()
    s.configuration.step = i
    s.particles.N = 5
    s.particles.types = ["A", "B"]
    s.particles.typeid = [0, 0, 1, 1, 1]
    s.particles.position = np.random.random(size=(5, 3))
    s.configuration.box = [3, 3, 3, 0, 0, 0]
    if bonds:
        s.bonds.N = 2
        s.bonds.types = ["AB"]
        s.bonds.typeid = [0, 0]
        s.bonds.group = [[0, 2], [1, 3]]
    if images:
        s.particles.image = np.full(shape=(5, 3), fill_value=i)
    else:
        s.particles.image = np.zeros(shape=(5, 3))
    s.validate()
    return s


def create_gsd(filename, bonds=False, images=False):
    with gsd.hoomd.open(name=filename, mode="wb") as f:
        f.extend(
            (create_frame(i, bonds=bonds, images=images) for i in range(10))
        )
