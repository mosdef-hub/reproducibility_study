import freud
import gsd
import gsd.hoomd
import numpy as np


def create_frame_random(i, bonds, n=50, L=10, seed=42):
    np.random.seed(seed)
    s = gsd.hoomd.Snapshot()
    s.configuration.step = i
    s.particles.N = n
    # First half of particles are "A" and the rest are "B"
    half = n - n // 2
    s.particles.types = ["A", "B"]
    s.particles.typeid = np.zeros(n)
    s.particles.typeid[half:] = 1
    s.particles.position = np.random.random(size=(n, 3)) * L - L / 2
    s.configuration.box = [L, L, L, 0, 0, 0]
    if bonds:
        s.bonds.types = ["AB"]
        # Create bonds between A-B particles within d_cut
        box = freud.box.Box.cube(L)
        aq = freud.locality.AABBQuery(box, s.particles.position[:half])
        d_cut = 2
        group = []
        selection = {"num_neighbors": 1}
        for i, j, d in aq.query(s.particles.position[half:], selection):
            if d < d_cut:
                group.append((j, i + half))
        s.bonds.typeid = np.zeros(len(group))
        s.bonds.group = group
        s.bonds.N = len(group)
    s.particles.image = np.zeros(shape=(n, 3))
    s.validate()
    return s


def create_frame_xstal(i, L=10, seed=42):
    box, points = freud.data.UnitCell.fcc().generate_system(
        num_replicas=L, sigma_noise=0.05, seed=seed
    )
    s = gsd.hoomd.Snapshot()
    s.configuration.step = i
    n = len(points)
    s.particles.N = n
    s.particles.types = ["A"]
    s.particles.typeid = np.zeros(n)
    s.particles.position = points
    s.configuration.box = [L, L, L, 0, 0, 0]
    s.particles.image = np.zeros(shape=(n, 3))
    s.validate()
    return s


def create_gsd(filename, bonds=False, system="random"):
    if system == "random":
        with gsd.hoomd.open(name=filename, mode="wb") as f:
            f.extend(
                [create_frame_random(i, bonds=bonds, seed=i) for i in range(10)]
            )
    else:
        with gsd.hoomd.open(name=filename, mode="wb") as f:
            f.extend([create_frame_xstal(i, seed=i) for i in range(10)])


if __name__ == "__main__":
    create_gsd("traj_random.gsd", bonds=False, system="random")
    create_gsd("traj_random_bonds.gsd", bonds=True, system="random")
    create_gsd("traj_xstal.gsd", system="xstal")
