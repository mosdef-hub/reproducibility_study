"""Methods used to create systems from job statepoint."""
import mbuild as mb
import signac
from mbuild import Box

from project.src import molecules

from .methane_ua import MethaneUA


def SystemBuilder(job):
    """Construct systems according to job statepoint."""
    # Update this dict as new reicpes are made
    molecule_dict = {
        "methaneUA": MethaneUA,
        "pentaneUA": None,
        "benzeneUA": None,
        "waterSPC/E": None,
        "ethanolAA": None,
    }
    molecule = molecule_dict[job.sp["molecule"]]
    vap_box = mb.Box()
    liq_box = mb.Box()

    vap_box_cpd = mb.fill_box(
        compound=molecule, n_compound=job.sp["n_compounds"], box=vap_box
    )
    liq_box_cpd = mb.fill_box(
        compound=molecule, n_compound=job.sp["n_compounds"], box=liq_box
    )

    system = mb.Compound([vap_box_cpd, vap_box_cpd])
    return system
