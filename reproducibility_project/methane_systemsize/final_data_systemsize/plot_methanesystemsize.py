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
# [N, MCCCS-MN	MCCCS-MN-SEM	MCCCS-MN_MOD	MCCCS-MN_MOD-SEM	LAMMPS	LAMMPS-SEM	LAMMPS_MOD	LAMMPS_MOD-SEM MCCCS-MN-T-1 MCCCS-MN-T-1-SEM]

data = {}
N = [450, 600, 900, 1800, 3600, 7200]

for i, n in enumerate(N):
    data[n] = raw_data[i, 1:]


fig, axs = plt.subplots(1, 1, sharex=False, sharey=False, figsize=(4.5, 4))
ax2 = axs

N = [450, 600, 900, 1800, 3600, 7200]
onebyN = [1000 / k for k in N]

Lammps = []
Lammps_err = []
MC3S = []
MC3S_err = []

MC3S_mod = []
MC3S_mod_err = []

MC3S_T_1 = []
MC3S_T_1_err = []

for n in N:
    Lammps.append(data[n][4])
    Lammps_err.append(1.96 * data[n][5])
    MC3S.append(data[n][0])
    MC3S_err.append(1.96 * data[n][1])
    MC3S_mod.append(data[n][2])
    MC3S_mod_err.append(1.96 * data[n][3])
    MC3S_T_1.append(data[n][6])
    MC3S_T_1_err.append(1.96 * data[n][7])

print("OK")

res = stats.linregress(onebyN, Lammps)
# print(f"R-squared: {res.rvalue**2:.6f}")
# ax2.plot(
#    onebyN,
#    res.intercept + res.slope * np.array(onebyN),
#    "--",
#    color=colors["LAMMPS"],
# )
ax2.errorbar(
    onebyN,
    1000 * np.array(Lammps),
    1000 * np.array(Lammps_err),
    capsize=4,
    label="LAMMPS",  # + " " + "$R^2$" + f"={res.rvalue**2:.2f}",
    marker=symbols["LAMMPS"],
    color=colors["LAMMPS"],
    markersize=3,
)


res = stats.linregress(onebyN, MC3S)
print(f"R-squared: {res.rvalue**2:.6f}")
# ax2.plot(
#    onebyN,
#    res.intercept + res.slope * np.array(onebyN),
#    "--",
#    color=colors["MCCCS-MN"],
# )
ax2.errorbar(
    onebyN,
    1000 * np.array(MC3S),
    1000 * np.array(MC3S_err),
    capsize=4,
    label="MCCCS-MN",  # + " " + "$R^2$" + f"={res.rvalue**2:.2f}",
    marker=symbols["MCCCS-MN"],
    color=colors["MCCCS-MN"],
    markersize=markersize,
)

res = stats.linregress(onebyN, MC3S_mod)
print(f"R-squared: {res.rvalue**2:.6f}")
# ax2.plot(
#    onebyN,
#    res.intercept + res.slope * np.array(onebyN),
#    "--",
#    color=colors["MCCCS-MN (MOD)"],
# )
ax2.errorbar(
    onebyN,
    1000 * np.array(MC3S_mod),
    1000 * np.array(MC3S_mod_err),
    capsize=4,
    label="MCCCS-MN(mod)",  # + " " + "$R^2$" + f"={res.rvalue**2:.2f}",
    marker=symbols["MCCCS-MN (MOD)"],
    color=colors["MCCCS-MN (MOD)"],
    markersize=markersize,
)
ax2.set_xlabel("$1000/N$", fontsize=xlabelfs)
ax2.set_ylabel(r"$\rho$, kg/m$^3$", fontsize=ylabelfs)


## T-1 scatter
# ax2.errorbar(1000/450, 1000*0.3763038544921875, 1000*2.0476215387896637e-05, label = "MCCCS-MN "+"$T_{-1}$", marker = "o", color= "magenta")
# ax2.errorbar(1000/900, 1000*0.3760736484375, 1000*1.4200439459421028e-05, marker = "o", color= "magenta")
# axes limits

ax2.errorbar(
    onebyN[0:5],
    1000 * np.array(MC3S_T_1)[0:5],
    1000 * np.array(MC3S_T_1_err)[0:5],
    capsize=4,
    label="MCCCS-MN $T_{\mathrm{-1}}$",  # + " " + "$R^2$" + f"={res.rvalue**2:.2f}",
    marker=symbols["MCCCS-MN-T-1"],
    color=colors["MCCCS-MN-T-1"],
    markersize=markersize,
)


ax2.tick_params(axis="y", labelsize=ytickfs)
ax2.tick_params(axis="x", labelsize=ytickfs)
props = dict(boxstyle="round", facecolor="none", alpha=1, ec="grey")
ax2.legend(frameon=True, ncol=1, fontsize=legendfs, labelspacing=0.05)


ax2.set_xlim([0.0001, 2.5])
ax2.set_ylim([375.6, 376.5])

ax2.yaxis.set_major_locator(plt.MaxNLocator(4))

# adjusting ticks

for ax in [ax2]:
    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))


plt.tight_layout()

fig.tight_layout()
plt.savefig("methane_systemsize.pdf", dpi=900)
