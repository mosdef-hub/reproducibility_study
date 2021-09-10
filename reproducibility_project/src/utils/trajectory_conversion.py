"""Convert Cassandra trajectories to gsd."""
import io
from pathlib import Path

import gsd.hoomd
import mbuild as mb
import numpy as np
import parmed


def cassandra2gsd(Hpath, xyzpath, gsdpath, species_list):
    """Convert Cassandra H and xyz files to a gsd trajectory file."""
    Hpath = Path(Hpath)
    xyzpath = Path(xyzpath)
    gsdpath = Path(gsdpath)

    nspecies = len(species_list)
    nmols_old = np.zeros(nspecies, dtype=int)

    with Hpath.open() as Hfile, xyzpath.open() as xyzfile, gsd.hoomd.open(
        str(gsdpath), "wb"
    ) as gsdfile:
        while Hfile.readline():
            with io.StringIO() as buff:
                for i in range(3):
                    buff.write(Hfile.readline())
                buff.seek(0)
                lmat = np.loadtxt(buff) * 0.1
            Hfile.readline()
            nspecies_in_box = int(Hfile.readline().strip())
            nmols = np.zeros(nspecies, dtype=int)
            for i in range(nspecies_in_box):
                mol_line_split = Hfile.readline().strip().split()
                nmols[int(mol_line_split[0]) - 1] = int(mol_line_split[1])
            natoms = int(xyzfile.readline().strip())
            step = int(xyzfile.readline().strip()[-1])
            with io.StringIO() as buff:
                for i in range(natoms):
                    buff.write(xyzfile.readline())
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
            types = list(set(types))
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
            gsdfile.append(s)
            nmols_old = nmols
