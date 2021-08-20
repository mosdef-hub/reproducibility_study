import mbuild as mb
import numpy as np

from reproducibility_project.src.molecules.methane_ua import MethaneUA
from reproducibility_project.tests.base_test import BaseTest


class TestMethaneUA(BaseTest):
    """Tests to ensure MethaneUA behaves as expected."""

    def test_creation(self):
        meth = MethaneUA()
        assert isinstance(meth, mb.Compound)
        assert np.all(np.isclose(meth.xyz, [[0.0, 0.0, 0.0]]))
        assert meth.name == "MethaneUA"
        assert [part.name == "_CH4" for part in meth.particles()]
        assert len(meth.labels) == 1
