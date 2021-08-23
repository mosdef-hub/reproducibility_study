import foyer
import mbuild as mb

from reproducibility_project.src.molecules.ethanol_aa import EthanolAA
from reproducibility_project.tests.base_test import BaseTest


class TestEthanolAA(BaseTest):
    """Test to ensure BenzeneUA behaves as expected."""

    def test_creation(self):
        etoh = EthanolAA()
        assert isinstance(etoh, mb.Compound)
        assert etoh.name == "EthanolAA"
        assert etoh.n_particles == 6
        assert etoh.n_bonds == 6
        assert len(etoh.labels) == 1

    def test_paramters(self):
        oplsaa = foyer.Forcefield(name="oplsaa")
        param_struct = oplsaa.apply(EthanolAA())

        expected_types = [
            "opls_135",
            "opls_140",
            "opls_154",
            "opls_155",
            "opls_157",
        ]
        for atom in param_struct.atoms:
            assert atom.atom_type.name in expected_types

        assert len(param_struct.bond_types) == 4
        assert len(param_struct.angle_types) == 5
        assert len(param_struct.rb_torsion_types) == 4
