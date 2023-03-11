"""Compare O-H bond-length distribution between LAMMPS, MCCCS-MN-flex, and MCCCS-MN."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

rc("font", **{"family": "serif", "serif": ["Times"]})
# rc('text', usetex=True)

# plot settings
ms = 8  # markersize
xtickfs = 11  # xtickfontsize
xlabelfs = 14  # xlabelfontsize
ylabelfs = 14  # ylabelfontsize
ytickfs = 11  # ytickfontsize
titlefs = 14  # title size
legendfs = 12
alpha = 0.2

mcccs_fix = "bl_analysis_data/mcccs_ethanolAA_NPT-fixOH_280.0K_101.325kPa_cutoff_hard_lrc_energy_pressure/O-H_hist.txt"
mcccs_flex = "bl_analysis_data/mcccs_ethanolAA_NPT-flexOH_280.0K_101.325kPa_cutoff_hard_lrc_energy_pressure/O-H_hist.txt"
lammps_fix = "bl_analysis_data/lammps-VU_ethanolAA_NPT-fixOH_280.0K_101.325kPa_cutoff_hard_lrc_energy_pressure/O-H_hist.txt"
lammps_flex = "bl_analysis_data/lammps-VU_ethanolAA_NPT-flexOH_280.0K_101.325kPa_cutoff_hard_lrc_energy_pressure/O-H_hist.txt"

mcccs_fix_data = np.genfromtxt(mcccs_fix)
mcccs_flex_data = np.genfromtxt(mcccs_flex)
lammps_fix_data = np.genfromtxt(lammps_fix)
lammps_flex_data = np.genfromtxt(lammps_flex)

plt.figure()

plt.plot(mcccs_flex_data[:, 0], mcccs_flex_data[:, 1], label="MCCCS-MN-flex")
plt.plot(lammps_flex_data[:, 0], lammps_flex_data[:, 1], label="LAMMPS")
plt.plot(lammps_fix_data[:, 0], lammps_fix_data[:, 1], label="LAMMPS-fixOH")
plt.ylim([0, 20])
plt.xlabel("Bond length" + r" ($\mathrm{\AA}$)", fontsize=xlabelfs)
plt.ylabel("Probability", fontsize=ylabelfs)
plt.xticks(fontsize=xtickfs)
plt.yticks(fontsize=ytickfs)
plt.tight_layout()
plt.fill_between(mcccs_flex_data[:, 0], mcccs_flex_data[:, 1], alpha=alpha)
plt.fill_between(lammps_flex_data[:, 0], lammps_flex_data[:, 1], alpha=alpha)
plt.fill_between(lammps_fix_data[:, 0], lammps_fix_data[:, 1], alpha=alpha)
plt.axvline(
    0.945,
    ls="--",
    color="r",
    label="FF bond length",
)


plt.legend(fontsize=legendfs)
plt.savefig("comp.png")
