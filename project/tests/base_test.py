import os

import foyer
import pytest

from project.src import xmls


class BaseTest:
    @pytest.fixture
    def spceff(self, name="spce.xml"):
        abs_path = os.path.dirname(os.path.abspath(xmls.__file__))
        return foyer.Forcefield(forcefield_files=str(abs_path) + name)
