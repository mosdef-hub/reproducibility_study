"""Methods used to create systems from job statepoint."""
import mbuild as mb
from methane_ua import MethaneUA
from pentane_ua import PentaneUA


def construct_system(sp):
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

    Returns
    -------
    system: mb.Compound
        The constructed mbuild Compound
    """
    system = mb.Compound()

    # Update this dict as new recipes are made
    molecule_dict = {
        "methaneUA": MethaneUA,
        "pentaneUA": PentaneUA,
        "benzeneUA": None,
        "waterSPC/E": None,
        "ethanolAA": None,
    }
    molecule = molecule_dict[sp["molecule"]]
    liq_box = mb.Box([sp["box_L_liq"]] * 3)
    filled_liq_box = mb.fill_box(
        compound=[molecule], n_compounds=[sp["N_liquid"]], box=liq_box
    )
    system.add(liq_box)

    if sp["box_L_vap"] and sp["N_vap"]:
        vap_box = mb.Box([sp["box_L_vap"]] * 3)
        filled_vap_box = mb.fill_box(
            compound=[molecule], n_compounds=[sp["N_vap"]], box=vap_box
        )

        # Translate the vapor box to be on top of the liquid box
        new_vap_box_center = filled_liq_box.center
        new_vap_box_center[2] += sp["box_L_vap"]
        system.add(filled_vap_box)

    return system
