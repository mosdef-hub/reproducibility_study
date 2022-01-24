import numpy as np
import pytest
from pymbar import testsystems
from pymbar.testsystems.timeseries import correlated_timeseries_example

from reproducibility_project.src.analysis.equilibration import (
    is_equilibrated,
    trim_non_equilibrated,
)
from reproducibility_project.tests.base_test import BaseTest


class TestEquilibration(BaseTest):
    def test_is_equilibrated(self, correlated_data_tau100_n10000):
        data = correlated_data_tau100_n10000
        assert not is_equilibrated(
            data, threshold_fraction=0.80, threshold_neff=100, strict=True
        )[0]
        assert not is_equilibrated(
            data, threshold_fraction=0.40, threshold_neff=100, strict=True
        )[0]
        assert is_equilibrated(
            data, threshold_fraction=0.10, threshold_neff=1, strict=True
        )[0]

        assert not is_equilibrated(
            correlated_data_tau100_n10000,
            threshold_fraction=0.10,
            threshold_neff=5000,
            strict=True,
        )[0]
        assert not is_equilibrated(
            correlated_data_tau100_n10000,
            threshold_fraction=0.10,
            threshold_neff=9999,
            strict=True,
        )[0]
        assert is_equilibrated(
            correlated_data_tau100_n10000,
            threshold_fraction=0.10,
            threshold_neff=10,
        )[0]

    def test_incorrect_threshold_fraction(self, correlated_data_tau100_n10000):
        with pytest.raises(
            ValueError, match=r"Passed \'threshold_fraction\' value"
        ):
            is_equilibrated(
                correlated_data_tau100_n10000, threshold_fraction=2.0
            )

        with pytest.raises(
            ValueError, match=r"Passed \'threshold_fraction\' value"
        ):
            is_equilibrated(
                correlated_data_tau100_n10000, threshold_fraction=-2.0
            )

    def test_incorrect_threshold_neff(self, correlated_data_tau100_n10000):
        data = correlated_data_tau100_n10000
        with pytest.raises(
            ValueError, match=r"Passed \'threshold_neff\' value"
        ):
            is_equilibrated(data, threshold_fraction=0.75, threshold_neff=0)
        with pytest.raises(
            ValueError, match=r"Passed \'threshold_neff\' value"
        ):
            is_equilibrated(data, threshold_fraction=0.75, threshold_neff=-1)

    def test_return_trimmed_data(self, correlated_data_tau100_n10000):
        data = correlated_data_tau100_n10000
        [new_a_t, t0, g, Neff] = trim_non_equilibrated(
            data, threshold_fraction=0.2, threshold_neff=10
        )
        assert np.shape(new_a_t)[0] < np.shape(data)[0]

    def test_trim_high_threshold(self, correlated_data_tau100_n10000):
        data = correlated_data_tau100_n10000
        with pytest.raises(ValueError, match=r"Data with a threshold_fraction"):
            [new_a_t, t0, g, Neff] = trim_non_equilibrated(
                data,
                threshold_fraction=0.98,
                strict=True,
            )
        with pytest.raises(
            ValueError,
            match=r"Data with a threshold\_fraction of 0\.75 and threshold\_neff 10000 is not equilibrated\!",
        ):
            [new_a_t, t0, g, Neff] = trim_non_equilibrated(
                data, threshold_fraction=0.75, threshold_neff=10000, strict=True
            )
