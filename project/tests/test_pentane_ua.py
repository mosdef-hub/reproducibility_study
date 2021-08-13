import mbuild as mb
import numpy as np
import pytest

from project.src.molecules.pentane_ua import PentaneUA
from project.tests.base_test import BaseTest


class TestPentaneUA(BaseTest):
    """Tests to ensure PentaneUA behaves as expected."""

    def test_creation(self):
        pent = PentaneUA()
        assert isinstance(pent, mb.Compound)
        assert pent.n_particles == 7
        assert pent.n_bonds == 6
        assert pent.name == "PentaneUA"
        assert all([p.name in ["_CH2", "_CH3"] for p in pent.particles()])
        assert len(pent.labels) == 1
