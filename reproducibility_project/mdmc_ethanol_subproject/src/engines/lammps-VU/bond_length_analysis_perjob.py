"""Setup for signac, signac-flow, signac-dashboard for running LAMMPS-VU simulations for the reproducibility study."""
import fileinput
import math
import os
import pathlib
import shutil
from glob import glob

import flow
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


class Rahman(DefaultSlurmEnvironment):
    """Subclass of DefaultPBSEnvironment for VU's Rahman cluster."""

    template = "rahman_lmp.sh"


@Project.operation.with_directives({"walltime": 200})
@Project.pre(lambda j: j.sp.engine == "lammps-VU")
@Project.pre(lambda j: j.isfile("trajectory-npt.gsd"))
@Project.post(lambda j: j.isfile("bl_C-C.txt"))
def save_bl_dist(job):
    """Save bl distribution for 4 types of bonds."""
    import os

    import mdtraj as md
    import numpy as np

    with job:
        if job.sp.molecule == "ethanolAA":
            filePath = "bl_C-C.txt"
            if os.path.exists(filePath):
                os.remove(filePath)
                print("{} deleted from {}".format(filePath, job))
            else:
                print("Can not delete the file as it doesn't exists")

            traj_filename = "trajectory-npt.gsd"
            traj = md.load(traj_filename, top="box.gro")
            atoms_per_molecule = 9
            total_atoms = 4500
            total_molecules = int(total_atoms / atoms_per_molecule)

            ## C-C bond length
            atom1 = "C"
            atom2 = "C"
            mean = 0.1529

            pairs = []

            for molecule_index in range(total_molecules):
                index1 = molecule_index * atoms_per_molecule
                index2 = molecule_index * atoms_per_molecule + 1
                # out = md.compute_distances(traj, np.array([[index1,index2]]))
                pairs.append([index1, index2])

            pairs = np.array(pairs)
            out = md.compute_distances(traj, pairs)
            bond_lengths = np.ravel(out) * 10
            np.savetxt("bl_C-C.txt", bond_lengths)

            atom1 = "C"
            atom2 = "O"
            mean = 0.141

            pairs = []

            for molecule_index in range(total_molecules):
                index1 = molecule_index * atoms_per_molecule + 1
                index2 = molecule_index * atoms_per_molecule + 2
                # out = md.compute_distances(traj, np.array([[index1,index2]]))
                pairs.append([index1, index2])

            pairs = np.array(pairs)
            out = md.compute_distances(traj, pairs)
            bond_lengths = np.ravel(out) * 10
            np.savetxt("bl_C-O.txt", bond_lengths)

            atom1 = "O"
            atom2 = "H"
            mean = 0.0945

            pairs = []

            for molecule_index in range(total_molecules):
                index1 = molecule_index * atoms_per_molecule + 2
                index2 = molecule_index * atoms_per_molecule + 3
                # out = md.compute_distances(traj, np.array([[index1,index2]]))
                pairs.append([index1, index2])

            pairs = np.array(pairs)
            out = md.compute_distances(traj, pairs)
            bond_lengths = np.ravel(out) * 10
            np.savetxt("bl_O-H.txt", bond_lengths)

            atom1 = "C"
            atom2 = "H"
            mean = 0.109

            pairs = []

            for molecule_index in range(total_molecules):
                indexC1 = molecule_index * atoms_per_molecule
                indexH7 = molecule_index * atoms_per_molecule + 6
                indexH8 = molecule_index * atoms_per_molecule + 7
                indexH9 = molecule_index * atoms_per_molecule + 8

                indexC2 = molecule_index * atoms_per_molecule + 1
                indexH5 = molecule_index * atoms_per_molecule + 4
                indexH6 = molecule_index * atoms_per_molecule + 5

                # out = md.compute_distances(traj, np.array([[index1,index2]]))
                pairs.append([indexC1, indexH7])
                pairs.append([indexC1, indexH8])
                pairs.append([indexC1, indexH9])
                pairs.append([indexC2, indexH5])
                pairs.append([indexC2, indexH6])

            pairs = np.array(pairs)
            out = md.compute_distances(traj, pairs)
            bond_lengths = np.ravel(out) * 10

            np.savetxt("bl_C-H.txt", bond_lengths)


if __name__ == "__main__":
    pr = Project()
    pr.main()
