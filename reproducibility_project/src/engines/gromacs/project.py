"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib

import flow
from flow import environments

from reproducibility_project.src.engine_input.gromacs import mdp


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


"""Setting progress label"""


@Project.label
def job_init(job):
    """Label for the initialization step."""
    return job.isfile("init.gro", "init.top", "nvt.mdp", "npt.mdp")


@Project.label
def em_grompp(job):
    """Label for the grompp_em step."""
    return job.isfile("em.tpr")


@Project.label
def em_completed(job):
    """Label for the gmx_em step."""
    return job.isfile("em.gro")


@Project.label
def nvt_grompp(job):
    """Label for the grompp_nvt step."""
    return job.isfile("nvt.tpr")


@Project.label
def nvt_completed(job):
    """Label for the gmx_nvt step."""
    return job.isfile("nvt.gro")


@Project.label
def npt_grompp(job):
    """Label for the grompp_npt step."""
    return job.isfile("npt.tpr")


@Project.label
def npt_completed(job):
    """Label for the gmx_npt step."""
    return job.isfile("npt.gro")


"""Setting up workflow operation"""


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "gromacs")
def init_job(job):
    """Initialize individual job workspace, including mdp and molecular init files."""
    from project.src.molecules.system_builder import construct_system

    with job:
        # Create a Compound and save to gro and top files
        system = construct_system(job.sp)
        system.save(filename="init.gro")
        system.save(filename="init.top", forcefield_name=job.sp.forcefield_name)

        # Modify mdp files according to job statepoint parameters
        mdp_abs_path = os.path.dirname(os.path.abspath(mdp.__file__))
        mdps = {
            "em": {
                "fname": "em.mdp",
                "template": f"{mdp_abs_path}/em_template.mdp.jinja",
                "data": dict(),
            },
            "nvt": {
                "fname": "nvt.mdp",
                "template": f"{mdp_abs_path}/nvt_template.mdp.jinja",
                "data": {"temp": job.sp.temperature},
            },
            "npt": {
                "fname": "npt.mdp",
                "template": f"{mdp_abs_path}/npt_template.mdp.jinja",
                "data": {
                    "temp": job.sp.temperature,
                    "refp": job.sp.pressure,
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
        return None


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "gromacs")
@Project.pre(lambda j: j.isfile("init.gro"))
@Project.pre(lambda j: j.isfile("init.top"))
@flow.cmd
def grompp_em(job):
    """Run GROMACS grompp for the energy minimization step."""
    with job:
        em_mdp_path = "em.mdp"
        msg = f"gmx grompp -f {em_mdp_path} -o em.tpr -c init.gro -p init.top --maxwarn 1"
        return msg


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "gromacs")
@Project.pre(lambda j: j.isfile("em.tpr"))
@flow.cmd
def gmx_em(job):
    """Run GROMACS mdrun for the energy minimization step."""
    with job:
        return _mdrun_str("em")


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "gromacs")
@Project.pre(lambda j: j.isfile("em.gro"))
@flow.cmd
def grompp_nvt(job):
    """Run GROMACS grompp for the nvt step."""
    with job:
        # nvt_mdp_path = "../../engine_input/gromacs/mdp/nvt.mdp"
        nvt_mdp_path = "nvt.mdp"
        msg = f"gmx grompp -f {nvt_mdp_path} -o nvt.tpr -c em.gro -p init.top --maxwarn 1"
        return msg


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "gromacs")
@Project.pre(lambda j: j.isfile("nvt.tpr"))
@flow.cmd
def gmx_nvt(job):
    """Run GROMACS mdrun for the nvt step."""
    with job:
        return _mdrun_str("nvt")


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "gromacs")
@Project.pre(lambda j: j.isfile("nvt.gro"))
@flow.cmd
def grompp_npt(job):
    """Run GROMACS grompp for the npt step."""
    with job:
        # npt_mdp_path = "../../engine_input/gromacs/mdp/npt.mdp"
        npt_mdp_path = "npt.mdp"
        msg = f"gmx grompp -f {npt_mdp_path} -o npt.tpr -c em.gro -p init.top --maxwarn 1"
        return msg


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt.tpr"))
@flow.cmd
def gmx_npt(job):
    """Run GROMACS mdrun for the npt step."""
    with job:
        return _mdrun_str("npt")


def _mdrun_str(op):
    """Output an mdrun string for arbitrary operation."""
    msg = f"gmx mdrun -v deffnm {op} -s {op}.tpr -cpi {op}.cpt "
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
    import jinja2
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
    breakpoint()
