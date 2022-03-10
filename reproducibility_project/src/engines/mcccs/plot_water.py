"""Plotting water results to debig why caluclated pressure is different from set pressure."""
import os

import matplotlib.pyplot as plt
import numpy as np


def read_trans_max(folder):
    """Read the maximum translation displacement from the config files."""
    os.chdir(folder)
    trans_max = []

    f = open("trans_max_disp.txt", "r")
    for line in f:
        trans_max.append(float(line))
    os.chdir("..")
    return trans_max


trans_max_280 = read_trans_max(
    "waterSPCE_NPT_280.0K_101.325kPa_cutoff_hard_lrc_energy_pressure"
)
trans_max_300 = read_trans_max(
    "waterSPCE_NPT_300.0K_101.325kPa_cutoff_hard_lrc_energy_pressure"
)
trans_max_320 = read_trans_max(
    "waterSPCE_NPT_320.0K_101.325kPa_cutoff_hard_lrc_energy_pressure"
)

x = range(1, 17)
fig, axs = plt.subplots(3, 1)
fig.suptitle("Translation max displacement", fontsize=16)

# plt.title("Translation max displacement")
plt.tick_params("x", labelbottom=False)
ax1 = axs[0]
# ax1 = plt.subplot(312)
ax1.scatter(x, trans_max_280)
ax1.tick_params("x", labelbottom=False)
ax1.set_title(r"$T$ = 280 K ")
ax2 = plt.subplot(312, sharex=ax1)
plt.scatter(x, trans_max_300)
plt.ylabel(r"Max displacement ($\mathrm{\AA}$)")

plt.tick_params("x", labelbottom=False)
ax2.set_title(r"$T$ = 300 K ")
plt.tick_params("x", labelbottom=False)
ax3 = plt.subplot(313, sharex=ax1)
ax3.set_title(r"$T$ = 320 K ")
plt.scatter(x, trans_max_320)
plt.tight_layout()
plt.xlabel("Seed")
plt.tight_layout()

plt.savefig("trans_max.png")
plt.close()


def read_rot_max(folder):
    """Read the maximum rotation displacement from the config files."""
    os.chdir(folder)
    rot_max = []

    f = open("rot_max_disp.txt", "r")
    for line in f:
        rot_max.append(float(line))
    os.chdir("..")
    return rot_max


rot_max_280 = read_rot_max(
    "waterSPCE_NPT_280.0K_101.325kPa_cutoff_hard_lrc_energy_pressure"
)
rot_max_300 = read_rot_max(
    "waterSPCE_NPT_300.0K_101.325kPa_cutoff_hard_lrc_energy_pressure"
)
rot_max_320 = read_rot_max(
    "waterSPCE_NPT_320.0K_101.325kPa_cutoff_hard_lrc_energy_pressure"
)

x = range(1, 17)
fig, axs = plt.subplots(3, 1)
fig.suptitle("Rotation max displacement", fontsize=16)

# plt.title("Translation max displacement")
ax1 = axs[0]
# ax1 = plt.subplot(312)
ax1.scatter(x, rot_max_280)
plt.tick_params("x", labelbottom=False)
ax1.set_title(r"$T$ = 280 K ")
ax2 = plt.subplot(312, sharex=ax1)
plt.scatter(x, rot_max_300)
plt.ylabel(r"Max displacement (radians)")
ax1.tick_params("x", labelbottom=False)

plt.tick_params("x", labelbottom=False)
ax2.set_title(r"$T$ = 300 K ")
ax3 = plt.subplot(313, sharex=ax1)
ax3.set_title(r"$T$ = 320 K ")
plt.scatter(x, rot_max_320)
plt.tight_layout()
plt.xlabel("Seed")
plt.tight_layout()

plt.savefig("rot_max.png")
plt.close()


def read_vol_max(folder):
    """Read the maximum volume displacement from the config files."""
    os.chdir(folder)
    vol_max = []

    f = open("vol_max_disp.txt", "r")
    for line in f:
        vol_max.append(float(line))
    os.chdir("..")
    return vol_max


vol_max_280 = read_vol_max(
    "waterSPCE_NPT_280.0K_101.325kPa_cutoff_hard_lrc_energy_pressure"
)
vol_max_300 = read_vol_max(
    "waterSPCE_NPT_300.0K_101.325kPa_cutoff_hard_lrc_energy_pressure"
)
vol_max_320 = read_vol_max(
    "waterSPCE_NPT_320.0K_101.325kPa_cutoff_hard_lrc_energy_pressure"
)

x = range(1, 17)
fig, axs = plt.subplots(3, 1)
fig.suptitle("Volume max displacement", fontsize=16)

# plt.title("Translation max displacement")
ax1 = axs[0]
# ax1 = plt.subplot(312)
ax1.scatter(x, vol_max_280)
plt.tick_params("x", labelbottom=False)
ax1.set_title(r"$T$ = 280 K ")
ax2 = plt.subplot(312, sharex=ax1)
plt.scatter(x, vol_max_300)
plt.ylabel(r"Max displacement ($\mathrm{\AA}^3$)")
ax1.tick_params("x", labelbottom=False)

plt.tick_params("x", labelbottom=False)
ax2.set_title(r"$T$ = 300 K ")
ax3 = plt.subplot(313, sharex=ax1)
ax3.set_title(r"$T$ = 320 K ")
plt.scatter(x, vol_max_320)
plt.tight_layout()
plt.xlabel("Seed")
plt.tight_layout()
plt.savefig("vol_max.png")
