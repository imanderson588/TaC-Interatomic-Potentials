import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import yaml

# --- Publication style settings ---
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "font.weight": "bold",
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.2,
    "lines.markersize": 4,
})

# --- Data loading ---
with open("gsf_111.yaml") as f:
    data_111 = yaml.load(f, Loader=yaml.FullLoader)

with open("gsf_110.yaml") as f:
    data_110 = yaml.load(f, Loader=yaml.FullLoader)

with open("gsf_112.yaml") as f:
    data_112 = yaml.load(f, Loader=yaml.FullLoader)

# Energy shift applied per atom during model training with FitSNAP
TA_E_SHIFT = 3.44766  # eV/atom
C_E_SHIFT = 1.24305  # eV/atom


# Conversion factor: eV to mJ/m^2
EV_TO_MJ_M2 = 16021.77


def calculate_E_shift(n_Ta, n_C):
    """Return total energy shift for a given number of Ta and C atoms."""
    return (n_Ta * TA_E_SHIFT) + (n_C * C_E_SHIFT)


# --- Process {111} ACE data ---
# 864.26 Å² = cross-sectional area of the {111} slip plane
AREA_111 = 69.1408
N_TA_111, N_C_111 = 240, 240

shift_111 = [d['shift'] for d in data_111]
raw_energy_111 = [
    d['pe'] - calculate_E_shift(N_TA_111, N_C_111) for d in data_111]
e0_111 = -raw_energy_111[0]
energy_111 = [((e + e0_111) / AREA_111) * EV_TO_MJ_M2 for e in raw_energy_111]

# --- Process {110} ACE data ---
# 706.43 Å² = cross-sectional area of the {110} slip plane
AREA_110 = 56.514519324
N_TA_110, N_C_110 = 80, 80

shift_110 = [d['shift'] for d in data_110]
raw_energy_110 = [
    d['pe'] - calculate_E_shift(N_TA_110, N_C_110) for d in data_110]
e0_110 = -raw_energy_110[0]
energy_110 = [((e + e0_110) / AREA_110) * EV_TO_MJ_M2 for e in raw_energy_110]

AREA_112 = 155.86
N_TA_112, N_C_112 = 0, 0

shift_112 = [d['shift']*np.sqrt(6) for d in data_112]
raw_energy_112 = [
    d['pe'] - calculate_E_shift(N_TA_112, N_C_112) for d in data_112]
e0_112 = -raw_energy_112[0]
energy_112 = [((e + e0_112) / AREA_112) * EV_TO_MJ_M2 for e in raw_energy_112]


# --- DFT reference data ---
# {111} DFT — 17.28 Å² per unit cell, energies in eV
AREA_DFT_111 = 17.28
displacement_dft_111 = [i * 0.05 for i in range(21)]
energy_dft_111_raw = [
    -525.06323373, -524.88473275, -524.40899587, -523.79122746, -523.19056648,
    -522.71930282, -522.39251957, -522.18534352, -522.07143977, -522.01026916,
    -521.99219536, -522.00982826, -522.07111012, -522.19129554, -522.39599035,
    -522.71999113, -523.19076522, -523.78641198, -524.40702389, -524.88500245,
    -525.06313516,
]
energy_dft_111 = [
    ((e + 525.06323373) / AREA_DFT_111) * 1.602e4 for e in energy_dft_111_raw
]

# {110} DFT — 14.12 Å² per unit cell, energies in eV
AREA_DFT_110 = 14.12
displacement_dft_110 = [i * 0.05 for i in range(21)]
energy_dft_110_raw = [
    -528.16, -528.01, -527.68, -527.32, -527.04,
    -526.87, -526.78, -526.75, -526.74439891, -526.74297734, -526.74217727, -526.74297734, -
    526.74439891, -526.75, -526.78, -526.87, -
    527.04, -527.32, -527.68, -528.01, -528.16
]
energy_dft_110 = [
    ((e + 528.16) / AREA_DFT_110) * 1.602e4 for e in energy_dft_110_raw
]

# --- Plotting ---
color_111 = '#0072B2'  # blue
color_110 = '#D55E00'  # orange

fig, ax = plt.subplots(figsize=(6, 3))

ax.plot(displacement_dft_111, energy_dft_111, "o-",
        color=color_111, label=r"DFT <110>$\{111\}$")
ax.plot(shift_111,            energy_111,     "s--",
        color=color_111, label=r"MEAM <110>$\{111\}$")
ax.plot(displacement_dft_110, energy_dft_110, "o-",
        color=color_110, label=r"DFT <110>$\{110\}$")
ax.plot(shift_110,            energy_110,     "s--",
        color=color_110, label=r"MEAM <110>$\{110\}$")

ax.plot(shift_112[0:16],            energy_112[0:16],     "s--",
        color=color_110, label=r"MEAM <112>$\{111\}$")


ax.set_xlabel("Fractional Displacement")
ax.set_ylabel(r"GSF Energy (mJ/m$^{-2}$)")
ax.legend(frameon=False)

fig.tight_layout()
fig.savefig("GSF_both.png", dpi=300)
plt.show()
