"""Utilities to load forcefields based on forcefield names."""
import foyer
import os


def load_ff(name: str = None,) -> foyer.Forcefield:
    """Based on a forcefield name, return a foyer.Forcefield object.

    For the reproducibility project, multiple forcefield types are expected based on the molecule of study at that statepoint.
    This will return a foyer.Forcefield object based on a naming convention defined in the init.py within the reproducibility_project.

    Parameters
    ----------
    name : str, default=None, optional
        Forcefield name to load.
    """

    if name in ["oplsaa", "trappe-ua"]:
        return foyer.Forcefield(name=name)
    elif name == "spce":
        from reproducibility_project.src import xmls

        ff_path = (
            str(os.path.dirname(os.path.abspath(xmls.__file__))) + "/spce.xml"
        )
        return foyer.Forcefield(forcefield_files=ff_path)
    elif name == "benzene-ua":
        from reproducibility_project.src import xmls

        ff_name = "benzene_trappe-ua_like.xml"
        ff_path = (
            str(os.path.dirname(os.path.abspath(xmls.__file__))) + "/" + ff_name
        )
        return foyer.Forcefield(forcefield_files=ff_path)
    else:
        raise ValueError(f"Unexpected forcefield name. Forcefield name {name} is not currently supported.")
