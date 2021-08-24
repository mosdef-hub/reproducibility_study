import pytest
from pymbar import testsystems

from reproducibility_project.src.analysis.equilibration import is_equilibrated
from reproducibility_project.src.analysis.sampler import _decorr_sampling
from reproducibility_project.tests.base_test import BaseTest


class TestSampler(BaseTest):
    def test_not_equilibrated(self):
        data = testsystems.correlated_timeseries_example(
            N=1000, tau=200, seed=432
        )
        with pytest.raises(ValueError):
            start, stop, step = _decorr_sampling(data, threshold=0.80)

    def test_equilibrated(self):
        data = testsystems.correlated_timeseries_example(
            N=1000, tau=200, seed=432
        )
        is_equil, prod_start, ineff = is_equilibrated(
            data, threshold=0.10, nskip=1
        )
        start, stop, step = _decorr_sampling(data, threshold=0.10)
        assert start >= prod_start
        assert step >= ineff
