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
def job_init(job):
    return (job.isfile('init.gro', 'init.top', 'nvt.mdp', 'npt.mdp'))

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
@Project.pre(lambda j: j.sp.simulation_engine == "gromcas")
def init_job(job):
    "Initialize individual job workspace, including mdp and molecular init files"
    import mbuild as mb
    from project.src.molecules.system_builder import SystemBuilder
    with job:
        # Create a Compound and save to gro and top files
        ff_path = '' # Fill in path to Trappe-UA xml
        system = SystemBuilder(job)
        system.save(filename='init.gro', combining_rule='geometric')
        system.save(filename='init.top', forcefield_name='trappe')

        # Modify mdp files according to job statepoint
        return None

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "gromcas")
@Project.pre(lambda j: j.isfile("init.gro"))
@Project.pre(lambda j: j.isfile("init.top"))
@flow.cmd
def grompp_em(job):
    with job:
        em_mdp_path = "../../engine_input/gromacs/mdp/em.mdp"
        msg = f"gmx grompp -f {em_mdp_path} -o em.tpr -c init.gro -p init.top --maxwarn 1"
        return msg

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "gromcas")
@Project.pre(lambda j: j.isfile("em.tpr"))
@flow.cmd
def gmx_em(job):
    with job:
        return _mdrun_str("em")

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "gromcas")
@Project.pre(lambda j: j.isfile("em.gro"))
@flow.cmd
def grompp_nvt(job):
    with job:
        #nvt_mdp_path = "../../engine_input/gromacs/mdp/nvt.mdp"
        nvt_mdp_path = "nvt.mdp"
        msg = f"gmx grompp -f {nvt_mdp_path} -o nvt.tpr -c em.gro -p init.top --maxwarn 1"
        return msg

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "gromcas")
@Project.pre(lambda j: j.isfile("nvt.tpr"))
@flow.cmd
def gmx_nvt(job):
    with job:
        return _mdrun_str("nvt")

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "gromcas")
@Project.pre(lambda j: j.isfile("em.gro"))
@flow.cmd
def grompp_npt(job):
    with job:
        #npt_mdp_path = "../../engine_input/gromacs/mdp/npt.mdp"
        npt_mdp_path = "npt.mdp"
        msg = f"gmx grompp -f {npt_mdp_path} -o npt.tpr -c em.gro -p init.top --maxwarn 1"
        return msg

@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "gromcas")
@Project.pre(lambda j: j.isfile("npt.tpr"))
@flow.cmd
def gmx_npt(job):
    with job:
        return _mdrun_str("npt")


def _mdrun_str(op):
    msg = f"gmx mdrun -v deffnm {op} -s {op}.tpr -cpi {op}.cpt "
    return msg

if __name__ == "__main__":
    pr = Project()
    pr.main()
    breakpoint()
