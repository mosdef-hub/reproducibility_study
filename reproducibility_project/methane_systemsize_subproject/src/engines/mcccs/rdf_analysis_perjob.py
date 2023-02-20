"""Setup for signac, signac-flow, signac-dashboard for running MCCCS-MN simulations for the reproducibility study."""
import fileinput
import math
import os
import pathlib
import shutil
from glob import glob

import flow
import freud
import mdtraj as md
from flow import FlowProject, environments
from flow.environment import DefaultSlurmEnvironment


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


class Metropolis(DefaultSlurmEnvironment):  # Grid(StandardEnvironment):
    """Subclass of DefaultSlurmEnvironment for Siepmann group cluster."""

    # metropolis.chem.umn.edu
    hostname_pattern = r".*\.chem\.umn\.edu"
    template = "metropolis.sh"


ex = Project.make_group(name="ex")


@Project.operation.with_directives({"walltime": 200})
@Project.pre(lambda j: j.sp.engine == "mcccs")
@Project.pre(lambda j: j.sp.ensemble != "NPT-exlarge")
def find_C_C_rdf(job):
    """Save rdf txt."""
    import os

    import mdtraj as md
    import numpy as np

    with job:
        if job.sp.molecule == "methaneUA":
            filePath = "C-C_rdf.txt"
            if os.path.exists(filePath):
                os.remove(filePath)
                print("{} deleted from {}".format(filePath, job))
            else:
                print("Can not delete the file as it doesn't exist")

            traj_filename = "trajectory-npt.gsd"
            traj = md.load(traj_filename, top="init1.mol2")
            traj = traj[::100]
            carbon_indices = traj.top.select("element C")
            bins = 250
            r_min = 0
            r_max = 2.5 / 2
            rdf_list = []
            freud_rdf = freud.density.RDF(bins=bins, r_min=r_min, r_max=r_max)
            for system in zip(
                np.asarray(traj.unitcell_vectors),
                traj.xyz[:, carbon_indices, :],
            ):
                freud_rdf.compute(system, reset=False)
            np.savetxt(
                "C-C_rdf.txt",
                np.vstack((freud_rdf.bin_centers, freud_rdf.rdf)).T,
            )
            np.savetxt(
                "C-C_cdf.txt",
                np.vstack((freud_rdf.bin_centers, freud_rdf.n_r)).T,
            )


if __name__ == "__main__":
    pr = Project()
    for job in pr.find_jobs():
        if job.sp.long_range_correction == None:
            pr.update_statepoint(
                job, {"long_range_correction": "None"}, overwrite=True
            )
    pr.main()
