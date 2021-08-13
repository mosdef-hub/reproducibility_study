import mbuild as mb
import numpy as np
import pytest

from project.src.analysis.rdf import gsd_rdf
from project.tests.base_test import BaseTest


class TestRDF(BaseTest):
    """Tests to ensure rdf behaves as expected."""

    def test_rdf(self, job_gsdfile):
        rdf = gsd_rdf(job_gsdfile)
