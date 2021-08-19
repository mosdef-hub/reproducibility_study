import freud
import mbuild as mb
import numpy as np
import pytest

from reproducibility_project.src.analysis.rdf import _gsd_rdf, gsd_rdf
from reproducibility_project.tests.base_test import BaseTest


class TestRDF(BaseTest):
    """Tests to ensure rdf behaves as expected."""

    def test_rdf(self, job_gsdfile):
        rdf = gsd_rdf(job_gsdfile)

    def test_gsdfile_random(self, gsdfile_random):
        rdf = _gsd_rdf(gsdfile_random)
        assert isinstance(rdf, freud.density.RDF)
        assert np.isclose(max(rdf.rdf), 1.5776997)
        assert len(rdf.rdf) == 50

    def test_gsdfile_xstal(self, gsdfile_xstal):
        rdf = _gsd_rdf(gsdfile_xstal)
        assert isinstance(rdf, freud.density.RDF)
        assert np.isclose(max(rdf.rdf), 2.5770662)
        assert len(rdf.rdf) == 50
