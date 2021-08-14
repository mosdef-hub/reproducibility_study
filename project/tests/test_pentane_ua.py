import foyer
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
        assert pent.n_particles == 5
        assert pent.n_bonds == 4
        assert pent.name == "PentaneUA"
        assert all([p.name in ["_CH2", "_CH3"] for p in pent.particles()])
        assert len(pent.labels) == 1

    def test_angles(self):
        pentane = PentaneUA()
        a = pentane.xyz[0]
        b = pentane.xyz[1]
        c = pentane.xyz[2]

        ba = a - b
        bc = c - b

        cosine_angle = np.dot(ba, bc) / (
            np.linalg.norm(ba) * np.linalg.norm(bc)
        )
        angle = np.arccos(cosine_angle)

        assert np.isclose(np.degrees(angle), 114.0)

    def test_lengths(self):
        pentane = PentaneUA()
        assert np.isclose(
            np.linalg.norm(pentane.xyz[3] - pentane.xyz[0]), 0.154
        )

    def test_parameters(self, trappe_ua):
        param_struct = trappe_ua.apply(PentaneUA())

        ch3 = trappe_ua.get_parameters("atoms", key="CH3_sp3")
        ch2 = trappe_ua.get_parameters("atoms", key="CH2_sp3")
        ch2_ch2_bond = trappe_ua.get_parameters(
            "harmonic_bonds", key=["CH2", "CH2"], keys_are_atom_classes=True
        )
        ch2_ch3_bond = trappe_ua.get_parameters(
            "harmonic_bonds", key=["CH3", "CH2"], keys_are_atom_classes=True
        )
        ch2_ch2_ch2_angle = trappe_ua.get_parameters(
            "harmonic_angles",
            key=["CH2", "CH2", "CH2"],
            keys_are_atom_classes=True,
        )
        ch3_ch2_ch2_angle = trappe_ua.get_parameters(
            "harmonic_angles",
            key=["CH3", "CH2", "CH2"],
            keys_are_atom_classes=True,
        )
        name_dict = {"CH3_sp3": 0, "CH2_sp3": 0}
        for atom in param_struct.atoms:
            name_dict[atom.atom_type.name] = name_dict[atom.atom_type.name] + 1

        assert name_dict["CH3_sp3"] == 2
        assert name_dict["CH2_sp3"] == 3

        assert "lorentz" == param_struct.combining_rule

        assert np.all(
            np.isclose(ch3["sigma"], param_struct.atoms[3].atom_type.sigma / 10)
        )
        assert np.all(
            np.isclose(
                ch3["epsilon"],
                param_struct.atoms[3].atom_type.epsilon * 4.184,
            )
        )
        assert np.all(np.isclose(ch3["charge"], param_struct.atoms[3].charge))

        assert np.all(
            np.isclose(ch2["sigma"], param_struct.atoms[0].atom_type.sigma / 10)
        )
        assert np.all(
            np.isclose(
                ch2["epsilon"],
                param_struct.atoms[0].atom_type.epsilon * 4.184,
            )
        )
        assert np.all(np.isclose(ch2["charge"], param_struct.atoms[0].charge))

        assert np.all(
            np.isclose(
                ch2_ch2_bond["k"],
                param_struct.bond_types[0].k * (4.184 * 2 / (0.1 ** 2)),
            )
        )
        assert np.all(
            np.isclose(
                ch2_ch2_bond["length"], param_struct.bond_types[0].req / 10
            )
        )

        assert np.all(
            np.isclose(
                ch2_ch2_ch2_angle["theta"],
                np.deg2rad(param_struct.angle_types[0].theteq).round(8),
            )
        )
        assert np.all(
            np.isclose(
                ch2_ch2_ch2_angle["k"],
                param_struct.angle_types[0].k * 4.184 * 2,
            )
        )
        assert np.all(
            np.isclose(
                ch3_ch2_ch2_angle["theta"],
                np.deg2rad(param_struct.angle_types[0].theteq).round(8),
            )
        )
        assert np.all(
            np.isclose(
                ch3_ch2_ch2_angle["k"],
                param_struct.angle_types[0].k * 4.184 * 2,
            )
        )
