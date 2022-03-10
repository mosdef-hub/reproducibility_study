"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib

import flow
import numpy as np
from flow import directives, environments

import reproducibility_project.templates.ndcrc


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


# ____________________________________________________________________________
"""Setting progress label"""


@Project.label
@Project.pre(lambda j: j.sp.engine == "cassandra")
def FinishedSPECalc(job):
    """Check if the log-spe.txt has been created."""
    return job.isfile("log-spe.txt")


# _____________________________________________________________________
"""Setting up workflow operation"""


def LoadSystemSnapShot(job):
    """Create initial configurations of the system statepoint."""
    import mbuild as mb

    pr = Project()
    snapshot_directory = (
        pathlib.Path(pr.root_directory()) / "src" / "system_snapshots"
    )
    molecule = job.sp.molecule
    molecule_filename = molecule + ".json"
    box = mb.load(str(snapshot_directory / molecule_filename))
    return box
    # __________________________________________________


@Project.operation
@Project.pre(lambda j: j.sp.engine == "cassandra")
@Project.post(FinishedSPCECalc)
@directives(omp_num_threads=4)
@flow.with_job
def CalculateEnergy(job):
    """Load onto a cluster and output the point energy for the snapshot."""
    import foyer
    import mosdef_cassandra as mc
    import pandas as pd
    import unyt as u

    molecule = job.sp.molecule

    compound = get_molecule(job.sp)
    ffname = job.sp.forcefield_name
    ff = load_ff(ffname)
    structure = ff.apply(compound)

    species_list = [structure]
    cutoff = job.sp.r_cut * u.nm

    box = LoadSystemSnapshot(job)

    system = mc.System([box], species_list, job.sp.N_liquid)
    moveset = mc.MoveSet("nvt", species_list)
    moveset.prob_rotate = 0.0
    moveset.prob_translate = 1.0
    nvtmoves.prob_regrow = 0.0
    nvtmoves.max_translate = [[0.0]]

    cass_cutoffs = {
        ("hard", "None"): "cut",
        ("hard", "energy_pressure"): "cut_tail",
        ("shift", "None"): "cut_shift",
    }
    cutoff_style = cass_cutoffs[
        (job.sp.cutoff_style, job.sp.long_range_correction)
    ]

    proplist = [
        "energy_total",
        "volume",
        "nmols",
        "pressure",
        "density",
        "mass_density",
        "energy_angle",
        "energy_dihedral",
        "energy_intravdw",
        "energy_intraq",
        "energy_inter",
        "energy_intervdw",
        "energy_lrc",
        "energy_interq",
        "energy_recip",
        "energy_self",
        "enthalpy",
    ]

    if any([abs(a.charge) > 0.0 for a in structure.atoms]):
        charge_style = "ewald"
    else:
        charge_style = "none"

    cutoff = job.sp.r_cut * u.nm

    mc.run(
        system=system,
        moveset=moveset,
        run_type="production",
        run_length=1,
        temperature=job.sp.temperature * u.K,
        properties=proplist,
        cutoff_style=cutoff_style,
        charge_style=charge_style,
        vdw_cutoff=cutoff,
        charge_cutoff=cutoff,
        run_name="nvt_spe",
        prop_freq=1,
        coord_freq=1,
        units="steps",
    )

    prpvec = np.loadtxt("nvt_spe.out.prp")[1:]
    prp = pd.Series(prpvec, index=proplist)
    spe_dict = {
        "total_energy": None,
        "potential_energy": prp.energy_total,
        "vdw_energy": prp.energy_intravdw
        + prp.energy_intervdw
        + prp.energy_lrc,
        "coul_energy": prp.energy_intraq
        + prp.energy_interq
        + prp.energy_self
        + prp.energy_recip,
        "pair_energy": prp.energy_intraq
        + prp.energy_interq
        + prp.energy_self
        + prp.energy_recip
        + prp.energy_intravdw
        + prp.energy_intervdw
        + prp.energy_lrc,
        "bonds_energy": None,
        "angles_energy": prp.energy_angle,
        "dihedrals_energy": prp.energy_dihedral,
        "tail_energy": prp.energy_lrc,
        "kspace_energy": prp.energy_recip,
    }
    spe_series = pd.Series(spe_dict)
    spe_log = pd.DataFrame(spe_series).T
    spe_log.to_csv("log-spe.txt", header=True, index=False, sep=",")

    # __________________________________________________


if __name__ == "__main__":
    pr = Project()
    pr.main()
