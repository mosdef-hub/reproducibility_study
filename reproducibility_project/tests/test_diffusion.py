import freud
import numpy as np

from reproducibility_project.src.analysis.diffusion import _gsd_msd, gsd_msd
from reproducibility_project.tests.base_test import BaseTest


class TestDiffusion(BaseTest):
    """Tests to ensure MSD/diffusion behaves as expected."""

    def test_diffusion(self, job_gsdfile):
        msd = gsd_msd(job_gsdfile)
        assert job_gsdfile.isfile("msd.txt")
        assert job_gsdfile.isfile("msd.png")
        assert job_gsdfile.doc.get("diffusion_coefficient") is not None

    def test_gsdfile_random(self, gsdfile_random):
        msd, timesteps = _gsd_msd(gsdfile_random, 0, 1, False)
        assert isinstance(msd, freud.msd.MSD)
        assert np.isclose(max(msd.msd), 50.0015816)
        assert len(msd.msd) == len(timesteps)
