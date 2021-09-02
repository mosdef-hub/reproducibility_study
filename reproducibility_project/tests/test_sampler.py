import pytest
from pymbar import testsystems

from reproducibility_project.src.analysis.equilibration import is_equilibrated
from reproducibility_project.src.analysis.sampler import _decorr_sampling
from reproducibility_project.tests.base_test import BaseTest


class TestSampler(BaseTest):
    def test_not_equilibrated(self, correlated_data_tau100_n10000):
        data = correlated_data_tau100_n10000
        with pytest.raises(ValueError):
            start, stop, step = _decorr_sampling(
                data, threshold_fraction=0.80, threshold_neff=100
            )

    def test_equilibrated(self, correlated_data_tau100_n10000):
        data = correlated_data_tau100_n10000
        is_equil, prod_start, ineff, Neff = is_equilibrated(
            data,
            threshold_fraction=0.10,
            threshold_neff=10,
            nskip=1,
        )
        start, stop, step, neff_samples = _decorr_sampling(
            data, threshold_fraction=0.10, threshold_neff=10
        )
        assert start >= prod_start
        assert step >= ineff
        assert neff_samples > 1
