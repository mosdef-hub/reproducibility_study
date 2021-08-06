"""Facilitate the calculation of the RDF from simulation data."""

import freud
import matplotlib.pyplot as plt

# Get system from simulation data
# this includes the points and the simulation box
# for a gsd file this can be as simple as:
#
# import gsd, gsd.hoomd
# with gsd.hoomd.open(gsdfile) as trajectory:
#    system = trajectory[-1]
#
# however the system can be represented as simply as (box, points):
#
# points = np.random.rand(50,3)
# L = 10 # coordinates go from -L/2 to L/2
# points = points * L - L/2 # Shift box points
# box = [L, L, L]
# system = (box, points)

# Decide on reasonable values
bins = 50
r_min = 0.5
r_max = 2.5

# compute RDF
rdf = freud.density.RDF(bins=bins, r_min=r_min, r_max=r_max)
rdf.compute(system, reset=False)

# plot
fig, ax = plt.subplots()
ax.plot(rdf.bin_centers, rdf.rdf)
ax.set_xlabel("$r$")
ax.set_ylabel("$g(r)$")
ax.set_title("RDF")
