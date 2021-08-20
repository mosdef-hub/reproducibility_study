"""Test file to ensure that forcefield loading is behaving as expected."""
import pytest
from reproducibility_project.tests.base_test import BaseTest
from foyer import Forcefield
from reproducibility_project.src.utils.forcefields import load_ff

class TestForcefields(BaseTest):
    def test_nonetype(self):
        with pytest.raises(ValueError, match=r"Unexpected forcefield name"):
            load_ff()

    def test_incorrect_value(self):
        with pytest.raises(ValueError, match=r"Unexpected forcefield name"):
            load_ff(name='foo')

    @pytest.mark.parametrize("ff_name", [
        ("oplsaa"), ("trappe-ua"), ("spce"), ("benzene-ua")])
    def test_correct_ff_names(self, ff_name):
        ff = load_ff(name=ff_name)
        assert isinstance(ff, Forcefield)
