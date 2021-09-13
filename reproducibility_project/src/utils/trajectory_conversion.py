"""Convert trajectory files to gsd format."""
import io
from pathlib import Path

import gsd.hoomd
import hoomd
import mbuild as mb
import numpy as np
import parmed


def cassandra2gsd(h_path, xyz_path, gsd_path, species_list):
    """Convert Cassandra H and xyz files to a gsd trajectory file.

    Inputs:
        h_path: path-like object (such as string or pathlib.Path) containing
            the path to the Cassandra .H file containing the box dimensions.
        xyz_path: path-like object (such as string or pathlib.Path) containing
            the path to the Cassandra .xyz file containing the trajectory atom coordinates.
        gsd_path: path-like object (such as string or pathlib.Path) containing
            the path to the gsd file to be written.
        species_list: list of parameterized single-molecule parmed Structure objects
            with one element per species.  This should be the same as the species_list
            supplied to the MoSDeF Cassandra System and MoveSet objects used to generate
            the trajectory.

    """
    h_path = Path(h_path)
    xyz_path = Path(xyz_path)
    gsd_path = Path(gsd_path)

    nspecies = len(species_list)
    nmols_old = np.zeros(nspecies, dtype=int)

    with h_path.open() as h_file, xyz_path.open() as xyz_file, gsd.hoomd.open(
        gsd_path, "wb"
    ) as gsd_file:
        while h_file.readline():
            with io.StringIO() as buff:
                for i in range(3):
                    buff.write(h_file.readline())
                buff.seek(0)
                lmat = np.loadtxt(buff) * 0.1
            h_file.readline()
            nspecies_in_box = int(h_file.readline().strip())
            nmols = np.zeros(nspecies, dtype=int)
            for i in range(nspecies_in_box):
                mol_line_split = h_file.readline().strip().split()
                nmols[int(mol_line_split[0]) - 1] = int(mol_line_split[1])
            natoms = int(xyz_file.readline().strip())
            step = int(xyz_file.readline().strip()[-1])
            with io.StringIO() as buff:
                for i in range(natoms):
                    buff.write(xyz_file.readline())
                buff.seek(0)
                xyz = np.loadtxt(buff, usecols=(1, 2, 3)) * 0.1
            if any(nmols != nmols_old):
                typed_system = parmed.Structure()
                for i, parspec in enumerate(species_list):
                    n = nmols[i]
                    if n > 0:
                        typed_system += parspec * n
            bonds = [
                (bond.atom1.idx, bond.atom2.idx) for bond in typed_system.bonds
            ]
            all_types = [a.type for a in typed_system.atoms]
            types = list(set(all_types))
            s = gsd.hoomd.Snapshot()
            s.configuration.step = step
            s.particles.N = natoms
            s.particles.position = xyz
            s.particles.types = types
            s.particles.typeid = [types.index(i) for i in all_types]
            s.bonds.N = len(bonds)
            s.bonds.group = bonds
            # must be upper triangular
            # todo: verify whether matrix needs to be transposed for non-ortho boxes
            box = hoomd.Box.from_matrix(lmat)
            s.configuration.box = [
                box.Lx,
                box.Ly,
                box.Lz,
                box.xy,
                box.xz,
                box.yz,
            ]
            s.validate()
            gsd_file.append(s)
            nmols_old = nmols
