"""
Track edge-dislocation position vs. timestep for every applied stress level
using OVITO's Dislocation Analysis (DXA), and plot all stress levels on one
figure to identify the Peierls (depinning) stress.

Expected dump layout (one directory per stress level, one file per snapshot):

    dump/dump_<stress>/dump.<timestep>

as produced in LAMMPS by:

    mkdir dump/dump_${stress}
    dump mydump all custom 10 dump/dump_${stress}/dump.* id type x y z

Outputs:
    disloc-position-vs-time.png   dislocation y-position vs. timestep,
                                  one curve per stress level (same axes)
    disloc-summary.png            net displacement and late-time velocity
                                  vs. applied stress

Notes:
  - Carbon (Type 2) atoms are deleted before DXA so interstitials don't
    interfere with the structural analysis.
  - Positions are unwrapped across the periodic y-boundary using a
    minimum-image convention between consecutive frames, so motion of many
    box lengths is captured correctly (valid as long as the dislocation
    moves less than half a box length between dumps).
  - If the dislocation is dissociated into Shockley partials, the reported
    position is the vertex-weighted mean over all segments, i.e. the center
    of the dissociated core.
"""

import glob
import os
import re

import numpy as np
import matplotlib.pyplot as plt

from ovito.io import import_file
from ovito.modifiers import (
    SelectTypeModifier,
    DeleteSelectedModifier,
    DislocationAnalysisModifier,
)

# ---------------------------------------------------------------------------
# Sampling stride (in LAMMPS timesteps).
# Only the first frame, every SAMPLE_EVERY-th timestep after it, and the last
# frame are analyzed and plotted. Set to None to use every dumped frame.
# Frames are discarded BEFORE the OVITO pipeline is built, so this also skips
# the (expensive) DXA computation on the unused frames.
SAMPLE_EVERY = 500
# ---------------------------------------------------------------------------


def build_pipeline(files):
    """Build an OVITO pipeline over an ordered list of dump files (frames)."""
    pipeline = import_file(files)

    # Remove interstitial carbon before structural analysis.
    pipeline.modifiers.append(SelectTypeModifier(
        operate_on="particles", property="Particle Type", types={2}))
    pipeline.modifiers.append(DeleteSelectedModifier())

    # Extract dislocation lines directly.
    pipeline.modifiers.append(DislocationAnalysisModifier(
        input_crystal_structure=DislocationAnalysisModifier.Lattice.FCC))

    return pipeline


def dislocation_y(data):
    """Vertex-weighted mean y-position of all dislocation line segments.

    Returns np.nan if DXA found no dislocation in this frame (can happen
    transiently during fast motion); nan frames are skipped when unwrapping.
    """
    ys, weights = [], []
    for seg in data.dislocations.segments:
        pts = np.asarray(seg.points)
        ys.append(pts[:, 1].mean())
        weights.append(len(pts))
    if not ys:
        return np.nan
    return float(np.average(ys, weights=weights))


def unwrap_series(y, box_length):
    """Unwrap a 1D position time series across a periodic boundary."""
    y = np.asarray(y, dtype=float)
    out = np.full_like(y, np.nan)
    last = None       # last unwrapped value
    last_raw = None   # last raw (wrapped) value
    for i, yi in enumerate(y):
        if np.isnan(yi):
            continue
        if last is None:
            out[i] = yi
        else:
            d = yi - last_raw
            d -= box_length * np.round(d / box_length)
            out[i] = last + d
        last, last_raw = out[i], yi
    return out


def collect_runs(root='dump'):
    """Map each stress level to its (timesteps, files), both sorted by timestep.

    Looks for directories  <root>/dump_<stress>/  containing files  dump.<timestep>
    """
    runs = {}
    for d in sorted(glob.glob(os.path.join(root, 'dump_*'))):
        if not os.path.isdir(d):
            continue
        m = re.search(r'dump_(\d+)$', os.path.basename(d))
        if not m:
            continue
        stress = int(m.group(1))

        entries = []
        for f in glob.glob(os.path.join(d, 'dump.*')):
            m2 = re.search(r'dump\.(\d+)$', os.path.basename(f))
            if m2:
                entries.append((int(m2.group(1)), f))
        entries.sort()

        if entries:
            timesteps, files = zip(*entries)
            runs[stress] = (list(timesteps), list(files))

    return dict(sorted(runs.items()))


def thin_by_timestep(timesteps, files, stride):
    """Keep the first frame, then every `stride` timesteps, plus the last frame.

    Selection is based on the actual timestep values in the filenames
    (relative to the first timestep of the run), so it works whether or not
    reset_timestep is used between stress levels, and regardless of dump
    frequency. Returns (timesteps, files) unchanged if stride is falsy.
    """
    if not stride:
        return timesteps, files
    t0, t_last = timesteps[0], timesteps[-1]
    keep_t, keep_f = [], []
    for t, f in zip(timesteps, files):
        if t == t0 or t == t_last or (t - t0) % stride == 0:
            keep_t.append(t)
            keep_f.append(f)
    return keep_t, keep_f


if __name__ == '__main__':
    runs = collect_runs()
    if not runs:
        raise SystemExit(
            "No dump files found matching dump/dump_<stress>/dump.<timestep>")

    trajectories = {}   # stress -> (timesteps, unwrapped y(t))

    for stress, (timesteps, files) in runs.items():
        timesteps, files = thin_by_timestep(timesteps, files, SAMPLE_EVERY)
        pipeline = build_pipeline(files)
        n_frames = pipeline.source.num_frames

        first = pipeline.compute(0)
        Ly = first.cell[1, 1]   # y box length (orthogonal cell assumed)

        ys = [dislocation_y(pipeline.compute(i)) for i in range(n_frames)]
        y_unwrapped = unwrap_series(ys, Ly)
        trajectories[stress] = (np.asarray(timesteps), y_unwrapped)

        valid = y_unwrapped[~np.isnan(y_unwrapped)]
        net = valid[-1] - valid[0] if len(valid) > 1 else 0.0
        print(f"stress {stress} MPa: {n_frames} frames "
              f"(timesteps {timesteps[0]}-{timesteps[-1]}), "
              f"net displacement {net:+.2f} A")

    # ---- Plot 1: average dislocation position vs. timestep, all stresses ----
    fig, ax = plt.subplots(figsize=(7.5, 5))
    cmap = plt.cm.viridis
    stresses = list(trajectories)
    for i, stress in enumerate(stresses):
        t, y = trajectories[stress]
        ax.plot(t, y,
                color=cmap(i / max(len(stresses) - 1, 1)),
                lw=1.5, label=f'{stress} MPa')
    ax.set_xlabel('Timestep')
    ax.set_ylabel(r'Average dislocation y-position ($\mathrm{\AA}$)')
    ax.set_title('Dislocation position vs. time for each applied stress')
    ax.legend(fontsize=8, ncol=2, title='Applied stress')
    fig.tight_layout()
    fig.savefig('disloc-position-vs-time.png', dpi=300)

    # ---- Plot 2: summary vs. stress -----------------------------------------
    net_disp, late_vel = [], []
    for stress in stresses:
        t, y = trajectories[stress]
        valid = ~np.isnan(y)
        tv, yv = t[valid], y[valid]
        if len(yv) < 2:
            net_disp.append(0.0)
            late_vel.append(0.0)
            continue
        net_disp.append(yv[-1] - yv[0])
        half = len(yv) // 2
        dt = tv[-1] - tv[half]
        late_vel.append((yv[-1] - yv[half]) / dt if dt > 0 else 0.0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 7), sharex=True)
    ax1.plot(stresses, net_disp, 'ko-')
    ax1.set_ylabel(r'Net displacement ($\mathrm{\AA}$)')
    ax2.plot(stresses, late_vel, 'ks-')
    ax2.set_ylabel(r'Late-time velocity ($\mathrm{\AA}$/timestep)')
    ax2.set_xlabel('Applied stress (MPa)')
    ax1.set_title('Depinning summary')
    fig.tight_layout()
    fig.savefig('disloc-summary.png', dpi=300)
    plt.show()
