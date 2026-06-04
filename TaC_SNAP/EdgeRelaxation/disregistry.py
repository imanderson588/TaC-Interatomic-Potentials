#!/usr/bin/env python3
"""
Compute and plot the disregistry of a dislocation from two LAMMPS data files.

    delta_i(xi) = u+_i(xi) - u-_i(xi)

u+ and u- are the per-atom displacements (defect minus base, matched by atom ID,
minimum-image on periodic axes) of the atomic planes just ABOVE and just BELOW
the slip plane, taken relative to the perfect-crystal base configuration.

Algorithm (after atomman):
  1. From the BASE config, identify atoms in the two planes neighboring the slip plane.
  2. For those atoms, compute u+ and u- (defect - base).
  3. Average any u+/u- sharing the same xi, then linearly interpolate both onto a
     common set of xi coordinates (the two planes need not share xi values).
  4. delta = u+ - u- on the common xi grid.

Only numpy + matplotlib required.
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use("Agg")

# ----------------------------------------------------------------------
# USER CONFIG  -- edit these for your geometry
# ----------------------------------------------------------------------
BASE_FILE = "perfect.lmp"      # perfect-crystal reference configuration
DEFECT_FILE = "EdgeRelaxed.lmp"  # dislocation configuration

NORMAL_AXIS = 2      # slip-plane normal:                 0=x, 1=y, 2=z
XI_AXIS = 1      # disregistry coordinate (glide dir): 0=x, 1=y, 2=z
SLIP_COORD = None   # slip-plane position along NORMAL_AXIS in Angstrom;
#   None -> midpoint of the model along the normal axis
PERIODIC = (True, True, False)  # periodicity per axis (x,y,z) for min-image
PLANE_TOL = 0.5    # +/- tolerance (Ang) for collecting atoms into a plane
XI_TOL = 0.05   # tolerance (Ang) for merging equal xi coordinates

COMPONENT_LABELS = ("x", "y", "z")
ATOM_TYPE = 1        # restrict disregistry to this atom type; None to use all types
OUT_PNG = "disregistry.png"
OUT_TXT = "disregistry.txt"
# ----------------------------------------------------------------------


def read_lammps_data(filename):
    """Minimal reader for an atom_style atomic LAMMPS data file.
    Returns ids[N], types[N], pos[N,3], box_lengths[3]."""
    with open(filename) as fh:
        lines = fh.readlines()

    natoms = None
    lo = [0.0, 0.0, 0.0]
    hi = [0.0, 0.0, 0.0]
    tilt = [0.0, 0.0, 0.0]
    atoms_at = None

    for i, raw in enumerate(lines):
        line = raw.split("#")[0].strip()
        if not line:
            continue
        tok = line.split()
        if len(tok) == 2 and tok[1] == "atoms":
            natoms = int(tok[0])
        elif "xlo" in line and "xhi" in line:
            lo[0], hi[0] = float(tok[0]), float(tok[1])
        elif "ylo" in line and "yhi" in line:
            lo[1], hi[1] = float(tok[0]), float(tok[1])
        elif "zlo" in line and "zhi" in line:
            lo[2], hi[2] = float(tok[0]), float(tok[1])
        elif "xy" in line and "xz" in line and "yz" in line:
            tilt = [float(tok[0]), float(tok[1]), float(tok[2])]
        elif tok[0] == "Atoms":
            atoms_at = i
            break

    if atoms_at is None or natoms is None:
        raise ValueError(
            f"Could not find Atoms section / atom count in {filename}")
    if any(abs(t) > 1e-8 for t in tilt):
        print(f"WARNING: {filename} is triclinic (tilt={tilt}); minimum-image "
              "uses orthogonal box lengths only -- verify periodic axes are orthogonal.")

    ids = np.empty(natoms, dtype=int)
    types = np.empty(natoms, dtype=int)
    pos = np.empty((natoms, 3), dtype=float)

    count = 0
    for raw in lines[atoms_at + 1:]:
        line = raw.split("#")[0].strip()
        if not line:
            if count == 0:
                continue          # blank line directly after "Atoms"
            break                 # blank line ends the Atoms block
        tok = line.split()
        ids[count] = int(tok[0])
        types[count] = int(tok[1])
        pos[count] = [float(tok[2]), float(tok[3]), float(tok[4])]
        count += 1
        if count == natoms:
            break

    if count != natoms:
        raise ValueError(
            f"Read {count} atoms, expected {natoms} in {filename}")

    box = np.array([hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]])
    return ids, types, pos, box


def find_neighbor_planes(coord_normal, slip_coord, tol):
    """Coordinates of the atomic planes immediately below and above slip_coord."""
    vals = np.sort(coord_normal)
    planes = [vals[0]]
    for v in vals[1:]:
        if v - planes[-1] > tol:
            planes.append(v)
    planes = np.array(planes)
    below = planes[planes < slip_coord]
    above = planes[planes > slip_coord]
    if below.size == 0 or above.size == 0:
        raise ValueError(
            f"Slip plane not bracketed by atomic planes; planes span "
            f"[{planes.min():.3f}, {planes.max():.3f}], slip_coord={slip_coord:.3f}")
    return below.max(), above.min()


def plane_profile(xi, disp, xi_tol):
    """Sort by xi and average displacements that share the same xi (within xi_tol)."""
    order = np.argsort(xi)
    xi, disp = xi[order], disp[order]
    grp_xi, grp_disp = [], []
    i, n = 0, len(xi)
    while i < n:
        j = i
        while j < n and xi[j] - xi[i] <= xi_tol:
            j += 1
        grp_xi.append(xi[i:j].mean())
        grp_disp.append(disp[i:j].mean(axis=0))
        i = j
    return np.array(grp_xi), np.array(grp_disp)


def main():
    b_ids, b_types, b_pos, b_box = read_lammps_data(BASE_FILE)
    d_ids, _, d_pos, _ = read_lammps_data(DEFECT_FILE)

    # match atoms one-to-one by atom ID
    b_sort, d_sort = np.argsort(b_ids), np.argsort(d_ids)
    if not np.array_equal(b_ids[b_sort], d_ids[d_sort]):
        raise ValueError("Atom IDs in the two files do not match one-to-one.")
    b_pos, d_pos = b_pos[b_sort], d_pos[d_sort]
    b_types = b_types[b_sort]

    # displacement, minimum-image on periodic axes
    du = d_pos - b_pos
    for ax in range(3):
        if PERIODIC[ax]:
            L = b_box[ax]
            du[:, ax] -= L * np.round(du[:, ax] / L)

    # restrict to the requested atom type
    if ATOM_TYPE is not None:
        mask = b_types == ATOM_TYPE
        print(f"Filtering to atom type {ATOM_TYPE}: {mask.sum()} of {len(mask)} atoms retained.")
        b_pos = b_pos[mask]
        du = du[mask]
    slip = SLIP_COORD
    if slip is None:
        slip = 0.5 * (b_pos[:, NORMAL_AXIS].min() +
                      b_pos[:, NORMAL_AXIS].max())

    below, above = find_neighbor_planes(b_pos[:, NORMAL_AXIS], slip, PLANE_TOL)
    print(f"slip plane (normal axis {NORMAL_AXIS}) at {slip:.3f} A;  "
          f"lower plane {below:.3f} A, upper plane {above:.3f} A")

    up_mask = np.abs(b_pos[:, NORMAL_AXIS] - above) < PLANE_TOL
    low_mask = np.abs(b_pos[:, NORMAL_AXIS] - below) < PLANE_TOL
    print(
        f"atoms in upper plane: {up_mask.sum()},  lower plane: {low_mask.sum()}")

    xi_up,  u_up = plane_profile(
        b_pos[up_mask,  XI_AXIS], du[up_mask],  XI_TOL)
    xi_low, u_low = plane_profile(
        b_pos[low_mask, XI_AXIS], du[low_mask], XI_TOL)

    # common xi grid: unique union restricted to the overlap of the two planes
    xi_common = np.unique(np.concatenate([xi_up, xi_low]))
    lo = max(xi_up.min(), xi_low.min())
    hi = min(xi_up.max(), xi_low.max())
    xi_common = xi_common[(xi_common >= lo) & (xi_common <= hi)]

    delta = np.empty((xi_common.size, 3))
    for c in range(3):
        u_p = np.interp(xi_common, xi_up,  u_up[:, c])
        u_m = np.interp(xi_common, xi_low, u_low[:, c])
        delta[:, c] = u_p - u_m

    header = "xi  " + \
        "  ".join(f"delta_{COMPONENT_LABELS[c]}" for c in range(3))
    np.savetxt(OUT_TXT, np.column_stack([xi_common, delta]), header=header)

    plt.figure(figsize=(7, 4.5))
    for c in range(3):
        plt.plot(xi_common, delta[:, c], "-o", ms=3,
                 label=fr"$\delta_{{{COMPONENT_LABELS[c]}}}$")
    plt.axhline(0, color="k", lw=0.5)
    plt.xlabel(fr"$\xi$ along {COMPONENT_LABELS[XI_AXIS]} [$\AA$]")
    plt.ylabel(r"disregistry $\delta_i$ [$\AA$]")
    plt.title("Dislocation disregistry")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=150)
    print(f"wrote {OUT_PNG} and {OUT_TXT}")


if __name__ == "__main__":
    main()
