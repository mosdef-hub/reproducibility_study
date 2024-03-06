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
import numpy as np
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


def calculate_force(r, sigma=3.73, epsilon=148, rcut=500):
    """Calculate force based on r."""
    if r >= rcut or r < 0.01:
        return 0
    term1 = (sigma**12) / (r**13)
    term2 = -0.5 * (sigma**6) / (r**7)
    return 48 * epsilon * (term1 + term2)


def calculate_jerk(r, sigma=3.73, epsilon=148, rcut=500):
    """Calculate jerk from r."""
    if r >= rcut or r < 0.01:
        return 0
    term1 = -13 * (sigma**12) / (r**14)
    term2 = 3.5 * (sigma**6) / (r**8)
    return 48 * epsilon * (term1 + term2)


def calculate_force_jerk_distance(distances):
    """Calculate force and jerk for a set of r."""
    forces = []
    jerks = []
    for distance in distances:
        forces.append(calculate_force(distance))
        jerks.append(calculate_jerk(distance))
    return np.sum(forces), np.sum(jerks)


@Project.operation.with_directives({"walltime": 200})
@Project.pre(lambda j: j.sp.engine == "mcccs")
def find_Tconfig(job):
    """Save Tconfig."""
    import os

    import mdtraj as md
    import numpy as np

    with job:
        if job.sp.molecule == "methaneUA":
            filePath = "Tconfig.txt"
            if os.path.exists(filePath):
                os.remove(filePath)
                print("{} deleted from {}".format(filePath, job))
            else:
                print("Can not delete the file as it doesn't exist")

            traj_filename = "trajectory-npt.gsd"
            traj = md.load(traj_filename, top="init1.mol2")
            trj = traj[::10]

            distance_matrix = np.zeros((trj.n_frames, trj.n_atoms, trj.n_atoms))
            force_sq_matrix = np.zeros((trj.n_frames, trj.n_atoms))
            jerk_matrix = np.zeros((trj.n_frames, trj.n_atoms))

            for i in range(trj.n_atoms):
                # print(i)
                # make list of lists
                pairs = []
                for j in range(trj.n_atoms):
                    pairs.append([i, j])

                r = 10 * md.compute_distances(trj, pairs, periodic=True)
                distance_matrix[:, i, :] = r

            for frame in range(trj.n_frames):
                # print(frame)
                for atom_no in range(trj.n_atoms):
                    distances = distance_matrix[frame, :, atom_no]
                    force, jerk = calculate_force_jerk_distance(distances)
                    # print(force, jerk)
                    force_sq_matrix[frame, atom_no] = force**2
                    jerk_matrix[frame, atom_no] = jerk

            # T_config = -1*np.average(force_sq_matrix)/np.average(jerk_matrix)
            T_config = -np.average(
                np.sum(force_sq_matrix, axis=1)
            ) / np.average(np.sum(jerk_matrix, axis=1))
            f = open("Tconfig.txt", "w")
            f.write(str(T_config))
            f.close()


if __name__ == "__main__":
    pr = Project()
    for job in pr.find_jobs():
        if job.sp.long_range_correction == None:
            pr.update_statepoint(
                job, {"long_range_correction": "None"}, overwrite=True
            )
    pr.main()
