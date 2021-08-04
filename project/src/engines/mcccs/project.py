"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib
import shutil
from glob import glob

import flow
from flow import environments


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


@Project.label
def has_fort_files(job):
    """Check if the job has all four equired fort.4 files."""
    return (
        job.isfile("fort.4.melt")
        and job.isfile("fort.4.cool")
        and job.isfile("fort.4.equil")
        and job.isfile("fort.4.prod")
    )


@Project.label
def has_restart_file(job):
    """Check if the job has a restart file."""
    return job.isfile("fort.77")


@Project.label
def has_restart_file(job):
    """Check if the job has a topmon (FF) file."""
    return job.isfile("fort.77")


@Project.label
def has_fort77maker(job):
    """Check if the job has a fort77maker file (obsolete)."""
    return os.path.isfile(
        Project().root_directory() + "/engines/mcccs/" + "fort77maker_onebox.py"
    )


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "mcccs")
@Project.post(has_fort_files)
def copy_files(job):
    """Copy the files for simulation from engine_input folder."""
    print(job.workspace())
    print(job.sp.molecule)
    for file in glob(
        Project().root_directory()
        + "/engine_input/mcccs/{}/fort.4.*".format(job.sp.molecule)
    ):
        shutil.copy(file, job.workspace() + "/")


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "mcccs")
@Project.post(has_fort77maker)
def copy_fort77maker(job):
    """Copy fort77maker_onebox.py from root directory to mcccs directory."""
    shutil.copy(
        Project().root_directory()
        + "/engine_input/mcccs/fort77maker_onebox.py",
        Project().root_directory() + "/engines/mcccs/",
    )


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "mcccs")
@Project.pre(has_fort77maker)
@Project.post(has_restart_file)
def make_restart_file(job):
    """Make a fort77 file for the job."""
    from fort77maker_onebox import fort77writer

    molecules = get_molecules(job)
    filled_box = get_system(job)
    fort77writer(molecules, filled_box, output_file="fort.77")


@Project.operation
@Project.pre(lambda j: j.sp.simulation_engine == "mcccs")
@Project.pre(has_restart_file)
@Project.pre(has_fort_files)
@project.pre(has_topmon)
def run_melt(job):
    """Run the melting stage."""
    return


if __name__ == "__main__":
    pr = Project()
    pr.main()
    print("The root dir is " + Project().root_directory())
    # breakpoint()
