"""Setup for signac, signac-flow, signac-dashboard for this study."""
# import foyer
import os
import pathlib

import flow
from flow import environments


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"

"""Setting progress label"""
@Project.label
def em_grompp(job):
    return job.isfile("em.tpr")

@Project.label
def em_completed(job):
    return job.isfile("em.gro")

@Project.label
def nvt_grompp(job):
    return job.isfile("nvt.tpr")

@Project.label
def nvt_completed(job):
    return job.isfile("nvt.gro")

@Project.label
def npt_grompp(job):
    return job.isfile("npt.tpr")

@Project.label
def npt_completed(job):
    return job.isfile("npt.gro")

"""Setting up workflow operation"""
@Project.operation
@flow.cmd
def grompp_em(job):
    _chdir(job)
    em_mdp_path = "../../engine_input/gromacs/mdp/em.mdp"
    msg = f"gmx grompp -f {em_mdp_path} -o em.tpr -c init.gro -p init.top --maxwarn 1"
    return msg

@Project.operation
@flow.cmd
def gmx_em(job):
    _chdir(job)
    return _mdrun_str("em")

@Project.operation
@flow.cmd
def grompp_nvt(job):
    _chdir(job)
    nvt_mdp_path = "../../engine_input/gromacs/mdp/nvt.mdp"
    msg = f"gmx grompp -f {nvt_mdp_path} -o nvt.tpr -c em.gro -p init.top --maxwarn 1"
    return msg

@Project.operation
@flow.cmd
def gmx_nvt(job):
    return _mdrun_str("nvt")

@Project.operation
@flow.cmd
def grompp_npt(job):
    _chdir(job)
    npt_mdp_path = "../../engine_input/gromacs/mdp/npt.mdp"
    msg = f"gmx grompp -f {npt_mdp_path} -o npt.tpr -c em.gro -p init.top --maxwarn 1"
    return msg

@Project.operation
@flow.cmd
def gmx_npt(job):
    _chdir(job)
    return _mdrun_str("npt")


def _mdrun_str(op):
    msg = f"gmx mdrun -v deffnm {op} -s {op}.tpr -cpi {op}.cpt "
    return msg

def _chdir(job):
    p = pathlib.Path(job.workspace())
    os.chdir(str(p.absolute()))

if __name__ == "__main__":
    pr = Project()
    pr.main()
    breakpoint()
