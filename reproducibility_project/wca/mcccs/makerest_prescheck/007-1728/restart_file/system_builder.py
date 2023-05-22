"""Methods used to create systems from job statepoint."""
import mbuild as mb
from mbuild.lib.molecules.water import WaterSPC

from reproducibility_project.src.molecules.benzene_ua import BenzeneUA
from reproducibility_project.src.molecules.ethanol_aa import EthanolAA
from reproducibility_project.src.molecules.methane_ua import MethaneUA
from reproducibility_project.src.molecules.pentane_ua import PentaneUA


def get_molecule(name):
    """Construct the mbuild molecule for the job statepoint.

    Parameters
    ----------
    sp: dict (from job.sp)
        Dictionary contains information necessary to construct a system.
        Stored as state of job. The dictionary should resemble:
        {"molecule": str,
         "engine": str,
         "replica": int,
         "temperature": float (in K),
         "pressure": float (in kPa),
         "ensemble": str,
         "N_liquid": int,
         "N_vap": int,
         "box_L_liq": int (nm),
         "box_L_vap", int (nm),
         "init_liq_den": float (g/cm3),
         "init_vap_den": float (g/cm3),
         "mass": float (g/mol),
         "forcefield_name": str,
         "cutoff_style": str,
         "r_cut": float (in nm)}

    Returns
    -------
    molecule
        Return mBuild molecule for the statepoint.
    """
    molecule_dict = {
        "methaneUA": MethaneUA(),
        "pentaneUA-flexible_bonds": PentaneUA(),
        "pentaneUA-constrain_bonds": PentaneUA(),
        "pentaneUA": PentaneUA(),
        "benzeneUA": BenzeneUA(),
        "waterSPCE": WaterSPC(),
        "ethanolAA": EthanolAA(),
    }
    molecule = molecule_dict["methaneUA"]
    molecule.name = "methaneUA"
    return molecule
