"""Setup for signac, signac-flow, signac-dashboard for running LAMMPS-VU simulations for the reproducibility study."""
import fileinput
import math
import os
import pathlib
import shutil
from glob import glob

import flow
import matplotlib.pyplot as plt
import mdtraj as md
from flow import FlowProject, environments
from flow.environment import DefaultSlurmEnvironment
from hbond_elliptical import calualateHBMap


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()
        self.data_dir = current_path.parents[1] / "data"
        self.ff_fn = self.data_dir / "forcefield.xml"


class Rahman(DefaultSlurmEnvironment):
    """Subclass of DefaultPBSEnvironment for VU's Rahman cluster."""

    template = "rahman_lmp.sh"


ex = Project.make_group(name="ex")


@Project.operation.with_directives({"walltime": 200})
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(lambda j: j.isfile("trajectory-npt.gsd"))
@Project.post(lambda j: j.isfile("n_hbond.txt"))
def hbond(job):
    """Save bl distribution for 4 types of bonds."""
    import os

    import mdtraj as md
    import numpy as np

    with job:
        if job.sp.molecule == "ethanolAA":
            filePath = "map_output.csv"
            traj_filename = "trajectory-npt.gsd"
            traj = md.load(traj_filename, top="box.gro")
            nbins_r = 200
            nbins_a = 200
            r_cutoff = 0.75
            skip_every_x_frames = 10
            sel_oxygen_head = "name O"
            sel_hydrogen = "name H"
            sel_oxygen_tail = "name O"
            list_names_hydrogen = ["H"]
            list_names_oxygen_head = ["O"]
            list_names_oxygen_tail = ["O"]
            (
                rdf_output,
                inter_output,
                map_output,
                hbond,
                hbond_time,
            ) = calualateHBMap(
                traj,
                r_cutoff,
                nbins_r,
                nbins_a,
                skip_every_x_frames,
                sel_oxygen_head,
                sel_oxygen_tail,
                sel_hydrogen,
                list_names_hydrogen,
                list_names_oxygen_head,
                list_names_oxygen_tail,
                bonded_pdb_provided=True,
            )
            inter_output[0] = inter_output[0] * 180 / np.pi
            plt.figure()
            cmap = plt.get_cmap("jet")
            plt.figure(figsize=(5, 3))
            plt.style.use("default")
            levels = np.linspace(0, 100, 11)
            cs = plt.contourf(
                rdf_output[0],
                inter_output[0],
                map_output,
                levels=levels,
                cmap=cmap,
            )
            plt.xlabel("r (nm)")
            plt.ylabel("\u03B8 (degrees)")
            plt.xlim([0.2, 0.4])
            plt.ylim([140, 180])
            plt.colorbar()

            plt.savefig("hbond.png")
            np.savetxt(
                os.path.join(job.ws, "map_output.csv"),
                map_output,
                delimiter=",",
            )
            np.savetxt("r.csv", rdf_output[0], delimiter=",")
            np.savetxt("theta.csv", inter_output[0], delimiter=",")
            print(np.mean(hbond_time))
            np.savetxt("n_hbond.txt", np.array([np.mean(hbond_time)]))


if __name__ == "__main__":
    pr = Project()
    pr.main()
