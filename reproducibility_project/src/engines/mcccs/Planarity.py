"""Compute benzene molecule planarity."""
import warnings

warnings.filterwarnings("ignore")
import matplotlib
import mbuild as mb
import mdtraj as md
import numpy as np
import parmed as pmd
import signac
from matplotlib import pyplot as plt

from reproducibility_project.src.molecules.system_builder import (
    construct_system,
)

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
}
def calculate_plane_of_best_fit(coords):
"""Compute plane of the benzene molecule."""
    # Calculate center of mass for uniform particle
    center_of_mass = np.mean(coords, axis=0)

    # Translate coordinates to center of mass
    coords -= center_of_mass

    # Calculate moment of inertia tensor
    inertia_tensor = np.zeros((3, 3))
    for i in range(coords.shape[0]):
        for j in range(3):
            for k in range(3):
                inertia_tensor[j][k] += coords[i][j] * coords[i][k]

    # Diagonalize moment of inertia tensor
    eigvals, eigvecs = np.linalg.eig(inertia_tensor)
    sort_indices = np.argsort(eigvals)

    # Select eigenvectors corresponding to smallest eigenvalues
    normal_vector = eigvecs[:, sort_indices[0]]

    # Normalize normal vector
    normal_vector /= np.linalg.norm(normal_vector)

    # Calculate plane coefficients
    a, b, c = normal_vector
    d = -np.dot(normal_vector, center_of_mass)

    return (a, b, c, d)

def distance_from_plane(point, plane_coeffs):
"""Compute particle distance from plane."""
    # Extract plane coefficients
    a, b, c, d = plane_coeffs

    # Calculate distance from point to plane
    numerator = abs(a*point[0] + b*point[1] + c*point[2] + d)
    denominator = np.sqrt(a**2 + b**2 + c**2)
    distance = numerator / denominator

    return distance

def rmsd_from_plane(points, plane_coeffs):
"""Compute rmsd of the particles from the plane."""
    # Calculate distances from points to plane
    distances = [distance_from_plane(point, plane_coeffs) for point in points]

    # Calculate sum of squared distances
    sum_squares = sum([distance**2 for distance in distances])

    # Calculate RMSD
    rmsd = np.sqrt(sum_squares / len(points))

    return rmsd


project = signac.get_project()
benzene_jobs = list()
for job in project.find_jobs({"molecule": "benzeneUA", "engine": "mcccs"}):
    print(job)
    benzene_jobs.append(job)

assert len(benzene_jobs) > 0
skip = 10
#  What are we plotting here
compiled_rmsd_npt = list()
compiled_rmsd_init = list()
for job in benzene_jobs:
    print(job)
    gsd_path = f"{job.ws}/trajectory-npt.gsd"
    mol2_path = f"{job.ws}/init.mol2"

    comp = construct_system(job.sp)[0]
    structure = comp.to_parmed()

    traj = md.load(gsd_path)
    traj = traj[::skip]
    print("This traj has {} frames".format(traj.n_frames))

    mapping = dict()
    for particle, atom in zip(
        comp.particles(), traj.topology.atoms):
        mapping[particle] = atom

    rmsds = list()
    original_rmsd = list()
    for child in comp.children:
        coords = [
            particle.pos for particle in child.particles()
        ]
        plane = calculate_plane_of_best_fit(coords)
        original_rmsd.append(rmsd_from_plane(coords, plane))
    compiled_rmsd_init.append(original_rmsd)

    for i in range(traj.n_frames):
        rmsds_per_frame = list()
        for child in comp.children:
            coords = [traj.xyz[i][mapping[particle].index]
                      for particle in child.particles()]
            plane = calculate_plane_of_best_fit(coords)
            rmsds_per_frame.append(rmsd_from_plane(coords, plane))
        rmsds.append(rmsds_per_frame)
    compiled_rmsd_npt.append(rmsds)

print(len(compiled_rmsd_npt))
print("##############################################")
print("The number of rmsds calculated is: ", len(compiled_rmsd_npt))
print("##############################################")
averaged_npt_rmsd = list()
yerr_npt_rmsd = list()
time_step = list()
for frame in range(traj.n_frames):
    frame_npt_rmsd = list()
    time_step.append(frame * 10*skip)
    print("Checking at frame number: ", frame)
    for i in range(16):
        frame_npt_rmsd.append(np.mean(compiled_rmsd_npt[i][frame]))
    averaged_npt_rmsd.append(np.mean(frame_npt_rmsd))
    yerr_npt_rmsd.append(np.std(frame_npt_rmsd))

averaged_init_rmsd = np.mean(compiled_rmsd_init)
yerr_init_rmsd = np.std(compiled_rmsd_init)

plt.plot(time_step,
         averaged_npt_rmsd,
        )

upper = np.array(averaged_npt_rmsd) + np.array(yerr_npt_rmsd)
lower = np.array(averaged_npt_rmsd) - np.array(yerr_npt_rmsd)
plt.fill_between(time_step,
                 lower,
                 upper,
                 alpha=0.2,
                 edgecolor='#1B2ACC',
                 facecolor='#089FFF',
                 label="Production",
                )

plt.axhline(averaged_init_rmsd,
            color = 'orange', linestyle = '--', label="Initial")

plt.xlim(0, 120000)
plt.ylabel("RMSD, $nm$", fontdict=font)
plt.xlabel("MCC", fontdict=font)
plt.legend(loc="center right")

plt.savefig("benzene_planarity/average_planarity.pdf", dpi=500, bbox_inches='tight')
