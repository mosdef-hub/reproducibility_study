"""Plot methane systemsize density."""
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("pdf")
import os
import shutil

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
from matplotlib import rc
from matplotlib.ticker import (
    FormatStrFormatter,
    MaxNLocator,
    MultipleLocator,
    NullFormatter,
    ScalarFormatter,
    StrMethodFormatter,
)
from scipy import stats

rc("font", **{"family": "serif", "serif": ["Times"]})
# rc('text', usetex=True)

# plot settings
ms = 8  # markersize
xtickfs = 11  # xtickfontsize
xlabelfs = 14  # xlabelfontsize
ylabelfs = 14  # ylabelfontsize
ytickfs = 11  # ytickfontsize
titlefs = 14  # title size
legendfs = 9
alpha = 0.2
markersize = 3

symbols = {}

symbols["MCCCS-MN"] = "^"
symbols["MCCCS-MN (MOD)"] = "<"

symbols["LAMMPS"] = "D"
symbols["MCCCS-MN-T-1"] = "o"

colors = {}
colors["MCCCS-MN"] = "b"
colors["MCCCS-MN (MOD)"] = "c"

colors["LAMMPS"] = "#7C1D6F"
colors["MCCCS-MN-T-1"] = "magenta"

engines = ["MCCCS-MN", "MCCCS-MN (MOD)", "LAMMPS", "MCCCS-MN-T-1"]

data_file = "methane_data.csv"

raw_data = np.genfromtxt(data_file, skip_header=1, delimiter=",")
# rcut, MCCCS-MN  ,    MCCCS-MN-CI,          LAMMPS,                LAMMPS-CI

data = {}
rcuts = [10, 14, 18, 24]

for i, rcut in enumerate(rcuts):
    data[rcut] = raw_data[i, 1:]

print(data)

fig, axs = plt.subplots(1, 1, sharex=False, sharey=False, figsize=(4.5, 4))
ax2 = axs


Lammps = []
Lammps_err = []
MC3S = []
MC3S_err = []


for rcut in rcuts:
    Lammps.append(data[rcut][2])
    Lammps_err.append(data[rcut][3])
    MC3S.append(data[rcut][0])
    MC3S_err.append(data[rcut][1])

print("OK")

ax2.errorbar(
    rcuts,
    1000 * np.array(Lammps),
    1000 * np.array(Lammps_err),
    capsize=4,
    label="MD",  # + " " + "$R^2$" + f"={res.rvalue**2:.2f}",
    marker=symbols["LAMMPS"],
    color=colors["LAMMPS"],
    markersize=3,
)


ax2.errorbar(
    rcuts,
    1000 * np.array(MC3S),
    1000 * np.array(MC3S_err),
    capsize=4,
    label="MC",  # + " " + "$R^2$" + f"={res.rvalue**2:.2f}",
    marker=symbols["MCCCS-MN"],
    color=colors["MCCCS-MN"],
    markersize=markersize,
)


ax2.tick_params(axis="y", labelsize=ytickfs)
ax2.tick_params(axis="x", labelsize=ytickfs)
props = dict(boxstyle="round", facecolor="none", alpha=1, ec="grey")
ax2.legend(frameon=True, ncol=1, fontsize=legendfs, labelspacing=0.05)


ax2.set_xlim([9, 25])
ax2.set_ylim([375, 376.5])

ax2.yaxis.set_major_locator(plt.MaxNLocator(4))

# adjusting ticks

for ax in [ax2]:
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))

plt.xticks([])

ax2.set_xlabel(
    "$r_{\mathrm{cut}}$" + ", " + "$\mathrm{\AA}$", fontsize=xlabelfs
)
ax2.set_ylabel(r"$\rho$, kg/m$^3$", fontsize=ylabelfs)
ax = plt.gca()
tcks = ["10", "14", "18", "24"]
locs = [10, 14, 18, 24]
plt.xticks(locs, tcks)
plt.tight_layout()
plt.minorticks_off()

fig.tight_layout()
plt.savefig("methane_rcut.pdf", dpi=900)
