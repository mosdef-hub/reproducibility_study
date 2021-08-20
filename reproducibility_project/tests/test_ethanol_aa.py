import foyer
import mbuild as mb

from reproducibility_project.src.molecules.ethanol_aa import EthanolAA
from reproducibility_project.tests.base_test import BaseTest


class TestEthanolAA(BaseTest):
    """Test to ensure BenzeneUA behaves as expected."""

    def test_creation(self):
        etoh = EthanolAA()
        assert isinstance(etoh, mb.Compound)
        assert etoh.name == "BenzeneUA"
        assert etoh.n_particles == 6
        assert etoh.n_bonds == 6
        assert len(etoh.labels) == 1

    def test_paramters(self):
        oplsaa = foyer.Forcefield(name="oplsaa")
        param_struct = oplsaa.apply(EthanolAA())
