"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib
import sys

import flow
import panedr
import unyt as u
from flow.environment import DefaultPBSEnvironment

from reproducibility_project.src.analysis.equilibration import is_equilibrated
from reproducibility_project.src.utils.forcefields import load_ff


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


class Rahman(DefaultPBSEnvironment):
    """Subclass of DefaultPBSEnvironment for VU's Rahman cluster."""

    template = "rahman_gmx.sh"

    @classmethod
    def add_args(cls, parser):
        """Add command line arguments to the submit call."""
        parser.add_argument(
            "--walltime",
            type=float,
            default=96,
            help="Walltime for this submission",
        )


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.post(lambda j: j.isfile("init.gro"))
@Project.post(lambda j: j.isfile("init.top"))
@Project.post(lambda j: j.isfile("em.mdp"))
@Project.post(lambda j: j.isfile("nvt.mdp"))
@Project.post(lambda j: j.isfile("npt_prod.mdp"))
@Project.post(lambda j: j.isfile("nvt_prod.mdp"))
@flow.with_job
def init_job(job):
    """Initialize individual job workspace, including mdp and molecular init files."""
    from reproducibility_project.src.engine_input.gromacs import mdp
    from reproducibility_project.src.molecules.system_builder import (
        construct_system,
    )

    # Create a Compound and save to gro and top files
    system = construct_system(job.sp)
    system[0].save(filename="init.gro", overwrite=True)
    ff = load_ff(job.sp.forcefield_name)
    param_system = ff.apply(system[0])
    param_system.save(
        "init.top",
        overwrite=True,
    )

    # Modify mdp files according to job statepoint parameters
    cutoff_styles = {"hard": "None", "shift": "Potential-shift"}

    pressure = job.sp.pressure * u.kPa
    mdp_abs_path = os.path.dirname(os.path.abspath(mdp.__file__))
    mdps = {
        "em": {
            "fname": "em.mdp",
            "template": f"{mdp_abs_path}/em_template.mdp.jinja",
            "data": {
                "r_cut": job.sp.r_cut,
                "cutoff_style": cutoff_styles[job.sp.cutoff_style],
                "temp": job.sp.temperature,
                "replica": job.sp.replica,
            },
        },
        "nvt": {
            "fname": "nvt.mdp",
            "template": f"{mdp_abs_path}/nvt_template.mdp.jinja",
            "data": {
                "nsteps": 2500000,
                "dt": 0.002,
                "temp": job.sp.temperature,
                "r_cut": job.sp.r_cut,
                "cutoff_style": cutoff_styles[job.sp.cutoff_style],
            },
        },
        "npt_prod": {
            "fname": "npt_prod.mdp",
            "template": f"{mdp_abs_path}/npt_template.mdp.jinja",
            "data": {
                "nsteps": 5000000,
                "dt": 0.001,
                "temp": job.sp.temperature,
                "refp": pressure.to_value("bar"),
                "r_cut": job.sp.r_cut,
                "cutoff_style": cutoff_styles[job.sp.cutoff_style],
            },
        },
        "nvt_prod": {
            "fname": "nvt_prod.mdp",
            "template": f"{mdp_abs_path}/nvt_template.mdp.jinja",
            "data": {
                "nsteps": 5000000,
                "dt": 0.001,
                "temp": job.sp.temperature,
                "r_cut": job.sp.r_cut,
                "cutoff_style": cutoff_styles[job.sp.cutoff_style],
            },
        },
    }

    for op, mdp in mdps.items():
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
@Project.post(lambda j: j.isfile("em.gro"))
@flow.with_job
@flow.cmd
def gmx_em(job):
    """Run GROMACS grompp for the energy minimization step."""
    em_mdp_path = "em.mdp"
    grompp = f"gmx_mpi grompp -f {em_mdp_path} -o em.tpr -c init.gro -p init.top --maxwarn 1"
    mdrun = _mdrun_str("em")
    return f"{grompp} && {mdrun}"


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("em.gro"))
@Project.post(lambda j: j.isfile("nvt.gro"))
@flow.with_job
@flow.cmd
def gmx_nvt(job):
    """Run GROMACS grompp for the nvt step."""
    nvt_mdp_path = "nvt.mdp"
    grompp = f"gmx_mpi grompp -f {nvt_mdp_path} -o nvt.tpr -c em.gro -p init.top --maxwarn 1"
    mdrun = _mdrun_str("nvt")
    return f"{grompp} && {mdrun}"


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("nvt.gro"))
@Project.post(lambda j: j.isfile("npt_prod.gro"))
@flow.with_job
@flow.cmd
def gmx_npt_prod(job):
    """Run GROMACS grompp for the npt step."""
    npt_mdp_path = "npt_prod.mdp"
    grompp = f"gmx_mpi grompp -f {npt_mdp_path} -o npt_prod.tpr -c nvt.gro -p init.top --maxwarn 1"
    mdrun = _mdrun_str("npt_prod")
    return f"{grompp} && {mdrun}"


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt_prod.gro"))
@Project.pre(lambda j: not equil_status(j, "npt_prod", "Potential"))
@Project.pre(lambda j: not equil_status(j, "npt_prod", "Volume"))
@Project.post(lambda j: equil_status(j, "npt_prod", "Potential"))
@Project.post(lambda j: equil_status(j, "npt_prod", "Volume"))
@flow.with_job
@flow.cmd
def extend_gmx_npt_prod(job):
    """Run GROMACS grompp for the npt step."""
    # Extend the npt run by 1000 ps (1 ns)
    extend = "gmx_mpi convert-tpr -s npt_prod.tpr -extend 1000 -o npt_prod.tpr"
    mdrun = _mdrun_str("npt_prod")
    return f"{extend} && {mdrun}"


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt_prod.gro"))
@Project.post(lambda j: j.isfile("nvt_prod.gro"))
@flow.with_job
@flow.cmd
def gmx_nvt_prod(job):
    """Run GROMACS grompp for the nvt step."""
    npt_mdp_path = "npt_prod.mdp"
    grompp = f"gmx_mpi grompp -f {npt_mdp_path} -o nvt_prod.tpr -c npt_prod.gro -p init.top --maxwarn 1"
    mdrun = _mdrun_str("nvt_prod")
    return f"{grompp} && {mdrun}"


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("nvt_prod.gro"))
@Project.pre(lambda j: not equil_status(j, "nvt_prod", "Potential"))
@Project.pre(lambda j: not equil_status(j, "nvt_prod", "Pressure"))
@Project.post(lambda j: equil_status(j, "nvt_prod", "Potential"))
@Project.post(lambda j: equil_status(j, "nvt_prod", "Pressure"))
@flow.with_job
@flow.cmd
def extend_gmx_nvt_prod_prod(job):
    """Run GROMACS grompp for the nvt step."""
    # Extend the npt run by 1000 ps (1 ns)
    extend = "gmx_mpi convert-tpr -s nvt_prod.tpr -extend 1000 -o nvt_prod.tpr"
    mdrun = _mdrun_str("nvt_prod")
    return f"{extend} && {mdrun}"


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt_prod.gro"))
@Project.pre(lambda j: equil_status(j, "npt_prod", "Potential"))
@Project.pre(lambda j: equil_status(j, "npt_prod", "Volume"))
@flow.with_job
def sample_npt_properties(job):
    """Sample properties of interest from npt edr."""
    import pandas as pd

    from reproducibility_project.src.analysis.sampler import (
        sample_job,
        write_subsampled_values,
    )

    p = pathlib.Path(job.workspace())
    data = panedr.edr_to_df(f"{str(p.absolute())}/npt.edr")
    # Properties of interest
    poi = {
        "Potential": "potential_energy",
        "Kinetic En.": "kinetic_energy",
        "Pressure": "pressure",
        "Temperature": "temperature",
        "Density": "density",
        "Volume": "volume",
    }

    tmp_df = pd.DataFrame()
    for idx, row in data.iterrows():
        tmp_df = tmp_df.append(row[list(poi.keys())])
    tmp_df.rename(poi, axis=1, inplace=True)
    tmp_df.insert(
        column="time_steps",
        value=[10000 * i for i in range(len(tmp_df))],
        loc=1,
    )
    tmp_df.to_csv("log-npt.txt", index=False, sep=" ")
    for prop in poi:
        sample_job(job, filename="log-npt.txt", variable=poi[prop])
        write_subsampled_values(
            job, property=poi[prop], property_filename="log-npt.txt"
        )


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("nvt_prod.gro"))
@Project.pre(lambda j: equil_status(j, "nvt_prod", "Potential"))
@Project.pre(lambda j: equil_status(j, "nvt_prod", "Pressure"))
@flow.with_job
def sample_nvt_properties(job):
    """Sample properties of interest from nvt edr."""
    import pandas as pd

    from reproducibility_project.src.analysis.sampler import (
        sample_job,
        write_subsampled_values,
    )

    p = pathlib.Path(job.workspace())
    data = panedr.edr_to_df(f"{str(p.absolute())}/nvt_prod.edr")
    # Properties of interest
    poi = {
        "Time": "time",
        "Potential": "potential_energy",
        "Kinetic En.": "kinetic_energy",
        "Pressure": "pressure",
        "Temperature": "temperature",
    }

    tmp_df = pd.DataFrame()
    for idx, row in data.iterrows():
        tmp_df = tmp_df.append(row[list(poi.keys())])
    tmp_df.rename(poi, axis=1, inplace=True)
    tmp_df.insert(
        column="time_steps",
        value=[10000 * i for i in range(len(tmp_df))],
        loc=1,
    )
    tmp_df.to_csv("log-nvt.txt", index=False, sep=" ")
    for prop in poi:
        sample_job(job, filename="log-nvt.txt", variable=poi[prop])
        write_subsampled_values(
            job, property=poi[prop], property_filename="log-nvt.txt"
        )


def _mdrun_str(op):
    """Output an mdrun string for arbitrary operation."""
    msg = f"gmx_mpi mdrun -v -deffnm {op} -s {op}.tpr -cpi {op}.cpt"
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


def equil_status(job, op, att):
    """Check equilibration status of specific attributes of specific operation."""
    p = pathlib.Path(job.workspace())
    if not job.isfile(f"{op}.edr"):
        return False
    else:
        data = panedr.edr_to_df(f"{str(p.absolute())}/{op}.edr")
        return is_equilibrated(data[att])[0]


if __name__ == "__main__":
    pr = Project()
    pr.main()
