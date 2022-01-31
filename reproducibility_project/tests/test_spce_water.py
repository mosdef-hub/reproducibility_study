import numpy as np
from mbuild.lib.molecules.water import WaterSPC

from reproducibility_project.tests.base_test import BaseTest


class TestSPCE(BaseTest):
    def test_distance(self):
        spce = WaterSPC()
        assert len([part for part in spce.particles()]) == 3
        parthw1 = [part for part in spce.particles_by_name("HW1")]
        parthw2 = [part for part in spce.particles_by_name("HW2")]
        partow1 = [part for part in spce.particles_by_name("OW")]
        assert np.all(
            np.isclose(np.linalg.norm(parthw1[0].xyz - partow1[0].xyz), 0.1)
        )
        assert np.all(
            np.isclose(np.linalg.norm(parthw2[0].xyz - partow1[0].xyz), 0.1)
        )
        assert np.all(
            np.isclose(
                np.linalg.norm(parthw1[0].xyz - parthw2[0].xyz).round(4),
                0.16330,
            )
        )

    def test_param_success(self, spceff):
        spceff.apply(WaterSPC())

    def test_parameters(self, spceff):
        param_struct = spceff.apply(WaterSPC())

        opls_116 = spceff.get_parameters("atoms", key="opls-116")
        opls_117 = spceff.get_parameters("atoms", key="opls-117")
        oh_bond = spceff.get_parameters(
            "harmonic_bonds", key=["HW", "OW"], keys_are_atom_classes=True
        )
        hoh_angle = spceff.get_parameters(
            "harmonic_angles",
            key=["HW", "OW", "HW"],
            keys_are_atom_classes=True,
        )

        assert "geometric" == param_struct.combining_rule
        assert np.all(
            np.isclose(
                opls_116["sigma"], param_struct.atoms[0].atom_type.sigma / 10
            )
        )
        assert np.all(
            np.isclose(
                opls_116["epsilon"],
                param_struct.atoms[0].atom_type.epsilon * 4.184,
            )
        )
        assert np.all(
            np.isclose(opls_116["charge"], param_struct.atoms[0].charge)
        )

        assert np.all(
            np.isclose(
                opls_117["sigma"], param_struct.atoms[1].atom_type.sigma / 10
            )
        )
        assert np.all(
            np.isclose(
                opls_117["epsilon"],
                param_struct.atoms[1].atom_type.epsilon * 4.184,
            )
        )
        assert np.all(
            np.isclose(opls_117["charge"], param_struct.atoms[1].charge)
        )

        assert np.all(
            np.isclose(
                opls_117["sigma"], param_struct.atoms[2].atom_type.sigma / 10
            )
        )
        assert np.all(
            np.isclose(
                opls_117["epsilon"],
                param_struct.atoms[2].atom_type.epsilon * 4.184,
            )
        )
        assert np.all(
            np.isclose(opls_117["charge"], param_struct.atoms[2].charge)
        )

        assert np.all(
            np.isclose(
                oh_bond["k"],
                param_struct.bond_types[0].k * (4.184 * 2 / (0.1**2)),
            )
        )
        assert np.all(
            np.isclose(oh_bond["length"], param_struct.bond_types[0].req / 10)
        )

        assert np.all(
            np.isclose(
                hoh_angle["theta"],
                np.deg2rad(param_struct.angle_types[0].theteq).round(8),
            )
        )
        assert np.all(
            np.isclose(
                hoh_angle["k"], param_struct.angle_types[0].k * 4.184 * 2
            )
        )
