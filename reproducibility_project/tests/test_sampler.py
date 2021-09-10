import pytest
from pymbar import testsystems

from reproducibility_project.src.analysis.equilibration import is_equilibrated
from reproducibility_project.src.analysis.sampler import (
    _decorr_sampling,
    write_subsampled_values,
)
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

    def test_write_subsampled_incorrect_params(self, tmp_job):
        with pytest.raises(
            TypeError,
            match=r"Expected input \'job\' of type signac\.contrib\.project\.Job",
        ):
            write_subsampled_values("foo", property="density", ensemble="npt")

        with pytest.raises(
            ValueError,
            match=r"Expected \'property\' to be a name of a property",
        ):
            write_subsampled_values(tmp_job, property="", ensemble="npt")

        with pytest.raises(
            ValueError,
            match=r"Attempting to overwrite already existing data for property",
        ):
            import numpy as np

            vals = [1, 2, 3, 4, 5, 6]
            out = f"""foo\n"""
            for val in vals:
                out = out + str(val) + "\n"
            with open(tmp_job.fn("log.txt"), "w") as fp:
                fp.write(out)
            tmp_job.data["npt/subsamples/foo"] = np.asarray([1, 2, 3, 4])
            tmp_job.doc["npt/sampling_results"] = {
                "foo": {"start": 1, "stop": 4, "step": 2, "Neff": 2}
            }

            write_subsampled_values(
                tmp_job, property="foo", ensemble="npt", overwrite=False
            )

    def test_file_missing(self, tmp_job):
        # by default, the tmp job is missing the file
        with pytest.raises(
            FileNotFoundError, match=r"File missing\.txt does not exist"
        ):
            tmp_job.doc["npt/sampling_results"] = {
                "foo": {"start": 1, "stop": 4, "step": 2, "Neff": 2}
            }
            write_subsampled_values(
                tmp_job,
                property="foo",
                ensemble="npt",
                property_filename="missing.txt",
                overwrite=False,
            )

    def test_correct_samples(self, tmp_job):
        import numpy as np

        vals = [1, 2, 3, 4, 5, 6]
        out = f"""foo\n"""
        for val in vals:
            out = out + str(val) + "\n"
        with open(tmp_job.fn("log.txt"), "w") as fp:
            fp.write(out)
        tmp_job.doc["npt/sampling_results"] = {
            "foo": {"start": 1, "stop": 4, "step": 1, "Neff": 4}
        }
        write_subsampled_values(
            tmp_job,
            property="foo",
            property_filename="log.txt",
            ensemble="npt",
            overwrite=False,
        )
        with tmp_job.data:
            assert len(tmp_job.data["npt/subsamples/foo"]) == 3
            np.testing.assert_array_equal(
                tmp_job.data["npt/subsamples/foo"], [2, 3, 4]
            )
