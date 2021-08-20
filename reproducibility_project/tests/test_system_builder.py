from reproducibility_project.src.molecules.system_builder import (
    construct_system,
)
from reproducibility_project.tests.base_test import BaseTest


class TestSystemBuilder(BaseTest):
    """Test to make sure the system_builder works as expected"""

    def test_liq_only(self, mock_job_npt):
        systems = construct_system(mock_job_npt)

    def test_liq_and_vap(self, mock_job_gemc):
        systems = construct_system(mock_job_gemc)
