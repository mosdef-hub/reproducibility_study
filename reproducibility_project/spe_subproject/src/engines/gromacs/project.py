"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib
import sys

import flow
import pandas as pd
import panedr
import unyt as u
from flow.environment import DefaultPBSEnvironment

from reproducibility_project.src.analysis.equilibration import is_equilibrated
from reproducibility_project.src.utils.forcefields import load_ff


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.post(lambda j: j.isfile("init.gro"))
@Project.post(lambda j: j.isfile("init.top"))
@Project.post(lambda j: j.isfile("nvt_ewald.mdp"))
@Project.post(lambda j: j.isfile("nvt_p3m.mdp"))
@flow.with_job
def LoadSystemSnapShot(job):
    """Create initial configurations of the system statepoint."""
    import mbuild as mb

    from reproducibility_project.spe_subproject.src.engine_input.gromacs import (
        mdp,
    )

    pr = Project()
    snapshot_directory = (
        pathlib.Path(pr.root_directory()) / "src" / "system_snapshots"
    )
    molecule = job.sp.molecule
    molecule_filename = molecule + ".json"
    box = mb.load(str(snapshot_directory / molecule_filename))
    box.save(filename="init.gro", precision=8, overwrite=True)
    # Apply forcefield and write out engine input files
    # __________________________________________________
    ff = load_ff(job.sp.forcefield_name)
    param_system = ff.apply(box)
    param_system.save(
        "init.top",
        overwrite=True,
    )

    # Modify mdp files according to job statepoint parameters
    cutoff_styles = {"hard": "None", "shift": "Potential-shift"}
    lrcs = {"None": "no", "energy_pressure": "EnerPres"}

    pressure = job.sp.pressure * u.kPa
    mdp_abs_path = os.path.dirname(os.path.abspath(mdp.__file__))
    mdps = {
        "nvt-ewald": {
            "fname": "nvt_ewald.mdp",
            "template": f"{mdp_abs_path}/nvt_template_ewald.mdp.jinja",
            "water-template": f"{mdp_abs_path}/nvt_template_water_ewald.mdp.jinja",
            "data": {
                "nsteps": 0,
                "dt": 0.001,
                "temp": job.sp.temperature,
                "refp": pressure.to_value("bar"),
                "r_cut": job.sp.r_cut,
                "cutoff_style": cutoff_styles[job.sp.cutoff_style],
                "lrc": lrcs[job.sp.long_range_correction],
            },
        },
        "nvt-p3m": {
            "fname": "nvt_p3m.mdp",
            "template": f"{mdp_abs_path}/nvt_template_p3m.mdp.jinja",
            "water-template": f"{mdp_abs_path}/nvt_template_water_p3m.mdp.jinja",
            "data": {
                "nsteps": 0,
                "dt": 0.001,
                "temp": job.sp.temperature,
                "r_cut": job.sp.r_cut,
                "cutoff_style": cutoff_styles[job.sp.cutoff_style],
                "lrc": lrcs[job.sp.long_range_correction],
            },
        },
    }

    for op, mdp in mdps.items():
        if job.sp.molecule == "waterSPCE":
            _setup_mdp(
                fname=mdp["fname"],
                template=mdp["water-template"],
                data=mdp["data"],
                overwrite=True,
            )
        else:
            _setup_mdp(
                fname=mdp["fname"],
                template=mdp["template"],
                data=mdp["data"],
                overwrite=True,
            )


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("init.gro"))
@Project.pre(lambda j: j.isfile("init.top"))
@Project.pre(lambda j: j.isfile("nvt_ewald.mdp"))
@Project.pre(lambda j: j.isfile("nvt_p3m.mdp"))
@Project.post(lambda j: j.isfile("nvt_ewald.edr"))
@Project.post(lambda j: j.isfile("nvt_p3m.edr"))
@flow.with_job
@flow.cmd
def CalculateEnergy(job):
    """Load onto a cluster and output the point energy for the snapshot."""
    nvt_ewald_mdp_path = "nvt_ewald.mdp"
    nvt_p3m_mdp_path = "nvt_p3m.mdp"
    grompp1 = f"gmx grompp -f {nvt_ewald_mdp_path} -o nvt_ewald.tpr -c init.gro -p init.top --maxwarn 1"
    mdrun1 = _mdrun_str("nvt_ewald")
    grompp2 = f"gmx grompp -f {nvt_p3m_mdp_path} -o nvt_p3m.tpr -c init.gro -p init.top --maxwarn 1"
    mdrun2 = _mdrun_str("nvt_p3m")

    return f"{grompp1} && {mdrun1} && {grompp2} && {mdrun2}"


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("nvt_ewald.edr"))
@Project.pre(lambda j: j.isfile("nvt_p3m.edr"))
@Project.post(lambda j: j.isfile("log-spe1.txt"))
@Project.post(lambda j: j.isfile("log-spe-ewald.txt"))
@Project.post(lambda j: j.isfile("log-spe-p3m.txt"))
@flow.with_job
def FormatTextFile(job):
    """Take the output from the simulation engine and convert it to log-spe.txt for data comparisons.

    See README.md for spe_subproject for formatting information.
    """
    import mdtraj
    import panedr

    p = pathlib.Path(job.workspace())
    data_ewald = panedr.edr_to_df(f"{str(p.absolute())}/nvt_ewald.edr").loc[0.0]
    data_p3m = panedr.edr_to_df(f"{str(p.absolute())}/nvt_p3m.edr").loc[0.0]

    to_drop = [
        "Vir-XX",
        "Vir-XY",
        "Vir-XZ",
        "Vir-YX",
        "Vir-YY",
        "Vir-YZ",
        "Vir-ZX",
        "Vir-ZY",
        "Vir-ZZ",
        "Pres-XX",
        "Pres-XY",
        "Pres-XZ",
        "Pres-YX",
        "Pres-YY",
        "Pres-YZ",
        "Pres-ZX",
        "Pres-ZY",
        "Pres-ZZ",
        "#Surf*SurfTen",
        "T-System",
        "Conserved En.",
        "Temperature",
        "Pres. DC",
        "Pressure",
        "Total Energy",
    ]
    for key in to_drop:
        data_ewald.pop(key)
        data_p3m.pop(key)

    spe_ewald = {
        "potential_energy": data_ewald.get("Potential", 0),
        "tot_vdw_energy": data_ewald.get("LJ (SR)", 0)
        + data_ewald.get("Disper. corr.", 0)
        + data_ewald.get("LJ-14", 0),
        "tail_energy": data_ewald.get("Disper. corr.", 0),
        "lj14": data_ewald.get("LJ-14", 0),
        "tot_electrostatics": data_ewald.get("Coulomb (SR)", 0)
        + data_ewald.get("Coul. recip.", 0)
        + data_ewald.get("Coulomb-14", 0),
        "short_range_electrostatics": data_ewald.get("Coulomb (SR)", 0),
        "long_range_electrostatics": data_ewald.get("Coul. recip."),
        "coulomb14": data_ewald.get("Coulomb-14", 0),
        "tot_pair_energy": data_ewald.get("LJ (SR)", 0)
        + data_ewald.get("Disper. corr.", 0)
        + data_ewald.get("LJ-14", 0)
        + data_ewald.get("Coulomb (SR)", 0)
        + data_ewald.get("Coul. recip.", 0)
        + data_ewald.get("Coulomb-14", 0),
        "bonds_energy": data_ewald.get("Bond", 0),
        "angles_energy": data_ewald.get("Angle", 0),
        "dihedrals_energy": data_ewald.get("Ryckaert-Bell.", 0)
        + data_ewald.get("Diehdral", 0),
        "tot_bonded_energy": data_ewald.get("Bond", 0)
        + data_ewald.get("Angle", 0)
        + data_ewald.get("Ryckaert-Bell.", 0)
        + data_ewald.get("Dihedral", 0),
        "intramolecular_energy": None,
        "intermolecular_energy": None,
    }
    spe_p3m = {
        "potential_energy": data_p3m.get("Potential", 0),
        "tot_vdw_energy": data_p3m.get("LJ (SR)", 0)
        + data_p3m.get("Disper. corr.", 0)
        + data_p3m.get("LJ-14", 0),
        "tail_energy": data_p3m.get("Disper. corr.", 0),
        "lj14": data_p3m.get("LJ-14", 0),
        "tot_electrostatics": data_p3m.get("Coulomb (SR)", 0)
        + data_p3m.get("Coul. recip.", 0)
        + data_p3m.get("Coulomb-14", 0),
        "short_range_electrostatics": data_p3m.get("Coulomb-14", 0),
        "long_range_electrostatics": data_p3m.get("Coul. recip."),
        "coulomb14": data_p3m.get("Coulomb-14", 0),
        "tot_pair_energy": data_p3m.get("LJ (SR)", 0)
        + data_p3m.get("Disper. corr.", 0)
        + data_p3m.get("LJ-14", 0)
        + data_p3m.get("Coulomb (SR)", 0)
        + data_p3m.get("Coul. recip.", 0)
        + data_p3m.get("Coulomb-14", 0),
        "bonds_energy": data_p3m.get("Bond", 0),
        "angles_energy": data_p3m.get("Angle", 0),
        "dihedrals_energy": data_p3m.get("Ryckaert-Bell.", 0)
        + data_p3m.get("Diehdral", 0),
        "tot_bonded_energy": data_p3m.get("Bond", 0)
        + data_p3m.get("Angle", 0)
        + data_p3m.get("Ryckaert-Bell.", 0)
        + data_p3m.get("Dihedral", 0),
        "intramolecular_energy": None,
        "intermolecular_energy": None,
    }

    spe_ewald_df = pd.DataFrame(spe_ewald, index=[0])
    spe_p3m_df = pd.DataFrame(spe_p3m, index=[0])
    spe_ewald_df.to_csv("log-spe-ewald.txt", header=True, index=False, sep=",")
    spe_p3m_df.to_csv("log-spe-p3m.txt", header=True, index=False, sep=",")

    spe_ewald_df.to_csv("log-spe.txt", header=True, index=False, sep=",")


"""
The below methods are adapted from
https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/drivers/gromacs.py
"""


def _get_gmx_energy_pair(gmx_energies):
    gmx_pairs = 0.0
    for key in ["LJ-14", "Coulomb-14"]:
        try:
            gmx_pairs += gmx_energies[key]
        except KeyError:
            pass
    return gmx_pairs


def _get_gmx_energy_torsion(gmx_energies):
    """Canonicalize torsion energies from a set of GROMACS energies."""
    gmx_torsion = 0.0
    for key in ["Torsion", "Ryckaert-Bell.", "Proper Dih."]:
        try:
            gmx_torsion += gmx_energies[key]
        except KeyError:
            pass

    return gmx_torsion


def _mdrun_str(op):
    """Output an mdrun string for arbitrary operation."""
    msg = f"gmx mdrun -v -deffnm {op} -s {op}.tpr -cpi {op}.cpt -nt 16"
    return msg


def _setup_mdp(fname, template, data, overwrite=False):
    """Create mdp files based on a template and provided data.

    Parameters
    ----------
    fname: str
        Name of the file to be saved out
    template: str, or jinja2.Template
        Either a jinja2.Template or path to a jinja template
    data: dict
        Dictionary storing data matched with the fields available in the template
    overwrite: bool, optional, default=False
        Options to overwrite (or not) existing mdp file of the

    Returns
    -------
    File saved with names defined by fname
    """
    from jinja2 import Template

    if isinstance(template, str):
        with open(template, "r") as f:
            template = Template(f.read())

    if not overwrite:
        if os.path.isfile(fname):
            raise FileExistsError(
                f"{fname} already exists. Set overwrite=True to write out."
            )

    rendered = template.render(data)
    with open(fname, "w") as f:
        f.write(rendered)

    return None


if __name__ == "__main__":
    pr = Project()
    pr.main()
