"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib
import sys

import flow


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@flow.with_job
def init_job(job):
    """Initialize individual job workspace, including mdp and molecular init files."""
    sys.path.append(Project().root_directory() + "/..")
    from reproducibility_project.src.engine_input.gromacs import mdp
    from reproducibility_project.src.molecules.system_builder import (
        construct_system,
    )

    # Create a Compound and save to gro and top files
    system = construct_system(job.sp)
    system[0].save(filename="init.gro", overwrite=True)
    if job.sp.forcefield_name in ["oplsaa", "trappe-ua"]:
        system[0].save(
            filename="init.top",
            forcefield_name=job.sp.forcefield_name,
            overwrite=True,
        )
    elif job.sp.forcefield_name == "spce":
        from reproducibility_project.src import xmls

        ff_path = (
            str(os.path.dirname(os.path.abspath(xmls.__file__))) + "/spce.xml"
        )
        system[0].save(
            filename="init.top", forcefield_files=ff_path, overwrite=True
        )
    elif job.sp.forcefield_name == "benzene-ua":
        from reproducibility_project.src import xmls

        ff_name = "benzene_trappe-ua_like.xml"
        ff_path = (
            str(os.path.dirname(os.path.abspath(xmls.__file__))) + "/" + ff_name
        )
        system[0].save(
            filename="init.top", forcefield_files=ff_path, overwrite=True
        )

    # Modify mdp files according to job statepoint parameters
    import unyt as u

    cutoff_styles = {"hard": "Cut-off"}
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
                "temp": job.sp.temperature,
                "r_cut": job.sp.r_cut,
                "cutoff_style": cutoff_styles[job.sp.cutoff_style],
            },
        },
        "npt": {
            "fname": "npt.mdp",
            "template": f"{mdp_abs_path}/npt_template.mdp.jinja",
            "data": {
                "temp": job.sp.temperature,
                "refp": pressure.to_value("bar"),
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
@flow.with_job
@flow.cmd
def grompp_em(job):
    """Run GROMACS grompp for the energy minimization step."""
    em_mdp_path = "em.mdp"
    msg = f"gmx grompp -f {em_mdp_path} -o em.tpr -c init.gro -p init.top --maxwarn 1"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("em.tpr"))
@flow.with_job
@flow.cmd
def gmx_em(job):
    """Run GROMACS mdrun for the energy minimization step."""
    return _mdrun_str("em")


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("em.gro"))
@flow.with_job
@flow.cmd
def grompp_nvt(job):
    """Run GROMACS grompp for the nvt step."""
    nvt_mdp_path = "nvt.mdp"
    msg = f"gmx grompp -f {nvt_mdp_path} -o nvt.tpr -c em.gro -p init.top --maxwarn 1"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("nvt.tpr"))
@flow.with_job
@flow.cmd
def gmx_nvt(job):
    """Run GROMACS mdrun for the nvt step."""
    return _mdrun_str("nvt")


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("nvt.gro"))
@flow.with_job
@flow.cmd
def grompp_npt(job):
    """Run GROMACS grompp for the npt step."""
    npt_mdp_path = "npt.mdp"
    msg = f"gmx grompp -f {npt_mdp_path} -o npt.tpr -c em.gro -p init.top --maxwarn 1"
    return msg


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt.tpr"))
@flow.with_job
@flow.cmd
def gmx_npt(job):
    """Run GROMACS mdrun for the npt step."""
    return _mdrun_str("npt")


def _mdrun_str(op):
    """Output an mdrun string for arbitrary operation."""
    msg = f"gmx mdrun -v -deffnm {op} -s {op}.tpr -cpi {op}.cpt "
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
