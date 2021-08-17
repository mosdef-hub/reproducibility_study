import numpy as np
import pytest
from pymbar import testsystems

from project.src.analysis.equlibration import (
    is_equilibrated,
    trim_non_equilibrated,
)

from .base_test import BaseTest


class TestEquilibration(BaseTest):
    def test_is_equilibrated(self):
        data = testsystems.correlated_timeseries_example(
            N=1000, tau=200, seed=432
        )
        assert not is_equilibrated(data, threshold=0.80)
        assert not is_equilibrated(data, threshold=0.40)
        assert is_equilibrated(data, threshold=0.10)

    def test_incorrect_threshold(self):
        data = testsystems.correlated_timeseries_example(
            N=1000, tau=200, seed=432
        )
        with pytest.raises(ValueError, match=r"Passed \'threshold\' value"):
            is_equilibrated(data, threshold=2.0)

        with pytest.raises(ValueError, match=r"Passed \'threshold\' value"):
            is_equilibrated(data, threshold=-2.0)

    def test_return_trimmed_data(self):
        data = testsystems.correlated_timeseries_example(
            N=1000, tau=200, seed=432
        )
        [new_a_t, g] = trim_non_equilibrated(data, threshold=0.2)
        assert np.shape(new_a_t)[0] < np.shape(data)[0]

    def test_trim_high_threshold(self):
        data = testsystems.correlated_timeseries_example(
            N=10000, tau=200, seed=432
        )
        with pytest.raises(ValueError, match=r"Data with a threshold"):
            [new_a_t, g, t0] = trim_non_equilibrated(data, threshold=0.98)
