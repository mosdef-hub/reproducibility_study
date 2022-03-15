import mbuild as mb
import numpy as np

from reproducibility_project.src.molecules.benzene_ua import BenzeneUA
from reproducibility_project.tests.base_test import BaseTest


class TestBenzeneUA(BaseTest):
    """Test to ensure BenzeneUA behaves as expected."""

    def test_creation(self):
        benz = BenzeneUA()
        assert isinstance(benz, mb.Compound)
        assert benz.name == "BenzeneUA"
        assert benz.n_particles == 6
        assert benz.n_bonds == 6
        assert [p.name == "_CH" for p in benz.particles()]
        assert len(benz.labels) == 1

    def test_parameters(self, benzene_ua_ff):
        param_struct = benzene_ua_ff.apply(BenzeneUA())

        nonbonded = benzene_ua_ff.get_parameters("atoms", key="CH_sp2")
        bond = benzene_ua_ff.get_parameters(
            "harmonic_bonds", key=["CH_E"] * 2, keys_are_atom_classes=True
        )
        angle = benzene_ua_ff.get_parameters(
            "harmonic_angles", key=["CH_E"] * 3, keys_are_atom_classes=True
        )
        rbtorsion = benzene_ua_ff.get_parameters(
            "rb_propers", key=["CH_E"] * 4, keys_are_atom_classes=True
        )

        for atom in param_struct.atoms:
            assert atom.atom_type.name == "CH_sp2"

        assert "lorentz" == param_struct.combining_rule
        assert np.isclose(
            nonbonded["sigma"], param_struct.atoms[0].atom_type.sigma / 10
        )
        assert np.isclose(
            nonbonded["epsilon"],
            param_struct.atoms[0].atom_type.epsilon * 4.184,
        )
        assert np.isclose(nonbonded["charge"], param_struct.atoms[0].charge)

        assert np.isclose(
            bond["k"], param_struct.bond_types[0].k * (4.184 * 2 / (0.1**2))
        )
        assert np.isclose(bond["length"], param_struct.bond_types[0].req / 10)

        assert np.isclose(
            angle["theta"],
            np.deg2rad(param_struct.angle_types[0].theteq).round(8),
        )
        assert np.isclose(angle["k"], param_struct.angle_types[0].k * 4.184 * 2)
