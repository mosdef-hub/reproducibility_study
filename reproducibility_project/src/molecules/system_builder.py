"""Methods used to create systems from job statepoint."""
import mbuild as mb
from mbuild.lib.molecules.water import WaterSPC

from reproducibility_project.src.molecules.benzene_ua import BenzeneUA
from reproducibility_project.src.molecules.ethanol_aa import EthanolAA
from reproducibility_project.src.molecules.methane_ua import MethaneUA
from reproducibility_project.src.molecules.pentane_ua import PentaneUA


def construct_system(sp, scale_liq_box=1.0, scale_vap_box=1.0, constrain=False):
    """Construct systems according to job statepoint.

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
    scale_liq_box: float, optional, default=1.0
        Option to scale sizes of the liquid box.
    scale_vap_box: float, optional, default=1.0
        Option to scale sizes of the vapor box.
    constrain: boolean, optional, default=False
        Option to use constrainmol on bond lengths to make them exactly equal to the bond lengths in the FF file.

    Returns
    -------
    [filled_liq_box, filled_vap_box]
        Return list of system as specified.
    """
    molecule = get_molecule(sp)
    liq_box = mb.Box([sp["box_L_liq"] * scale_liq_box] * 3)

    filled_liq_box = mb.fill_box(
        compound=[molecule], n_compounds=[sp["N_liquid"]], box=liq_box
    )

    if sp["box_L_vap"] and sp["N_vap"]:
        vap_box = mb.Box([sp["box_L_vap"] * scale_vap_box] * 3)
        filled_vap_box = mb.fill_box(
            compound=[molecule], n_compounds=[sp["N_vap"]], box=vap_box
        )
        boxes = [filled_liq_box, filled_vap_box]
    else:
        boxes = [filled_liq_box, None]

    if not constrain or sp.molecule == "methaneUA":
        return boxes

    # If we reached this far, we need constrainmol, foyer, and load_ff
    import foyer
    from constrainmol import ConstrainedMolecule

    from reproducibility_project.src.utils.forcefields import load_ff

    ff = load_ff(sp.forcefield_name)
    parmed_molecule = molecule.to_parmed()
    typed_molecule = ff.apply(parmed_molecule)
    constrain_mol = ConstrainedMolecule(typed_molecule)

    for box in boxes:
        if box is None:
            continue
        else:
            for mol in box.children:
                constrain_mol.update_xyz(mol.xyz * 10)  # nm to angstrom
                constrain_mol.solve()
                mol.xyz = constrain_mol.xyz / 10.0  # angstrom to nm

    return boxes


def get_molecule(sp):
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
        "pentaneUA": PentaneUA(),
        "benzeneUA": BenzeneUA(),
        "waterSPCE": WaterSPC(),
        "ethanolAA": EthanolAA(),
    }
    molecule = molecule_dict[sp["molecule"]]
    molecule.name = sp["molecule"]
    return molecule
