import foyer
import pytest


class BaseTest:
    @pytest.fixture
    def spceff(self, name="spce.xml"):
        return foyer.Forcefield(forcefield_files="../src/xmls/" + name)
