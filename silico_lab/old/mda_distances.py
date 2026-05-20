#!/usr/bin/env python3
"""Calculate and plot distances between selected atom pairs from a GROMACS trajectory."""

import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ── Atom pair definitions ────────────────────────────────────────────
# Each entry: (selection1, selection2, label)
PAIRS = [

    # loop C
    # 1-2 F200 ↔ T202 h-bond cap
    ("chainID A and resid 205 and name OH",  "chainID A and resid 202 and name OG1",
     "A:Tyr205-OH ↔ A:Thr202-OG1"),
    ("chainID C and resid 205 and name OH",  "chainID C and resid 202 and name OG1",
     "C:Tyr205-OH ↔ C:Thr202-OG1"),
    # 3-4 F200 ↔ F46 BS opening
    ("chainID A and resid 200 and name CG",  "chainID B and resid 46 and name CG",
     "A:200-CG ↔ B:46-CG"),
    ("chainID C and resid 200 and name CG",  "chainID D and resid 46 and name CG",
     "C:200-CG ↔ D:46-CG"),

    # GABA electro
    # 5-6 T202 ↔ ABU-C4
    ("chainID A and resid 202 and name OG1", "chainID H and resname ABU and name C4",
     "A:Thr202-OG1 ↔ H:ABU-C4"),
    ("chainID C and resid 202 and name OG1", "chainID I and resname ABU and name C4",
     "C:Thr202-OG1 ↔ I:ABU-C4"),
    # 7-8 R67 ↔ ABU-C4
    ("chainID B and resid 67 and name CZ",   "chainID H and resname ABU and name C4",
     "B:Arg67-CZ ↔ H:ABU-C4"),
    ("chainID D and resid 67 and name CZ",   "chainID I and resname ABU and name C4",
     "D:Arg67-CZ ↔ I:ABU-C4"),
    # 9-10 T130 ↔ ABU-C4
    ("chainID B and resid 130 and name OG1", "chainID H and resname ABU and name C4",
     "B:Thr130-OG1 ↔ H:ABU-C4"),
    ("chainID D and resid 130 and name OG1", "chainID I and resname ABU and name C4",
     "D:Thr130-OG1 ↔ I:ABU-C4"),
    # 11-12 E155 ↔ ABU-N
    ("chainID A and resid 155 and name CD",  "chainID H and resname ABU and name N",
     "A:Glu155-CD ↔ H:ABU-N"),
    ("chainID C and resid 155 and name CD",  "chainID I and resname ABU and name N",
     "C:Glu155-CD ↔ I:ABU-N"),
    # 13-14 Y97 ↔ ABU-N
    ("chainID A and resid 97 and name OH",   "chainID H and resname ABU and name N",
     "A:Tyr97-OH ↔ H:ABU-N"),
    ("chainID C and resid 97 and name OH",   "chainID I and resname ABU and name N",
     "C:Tyr97-OH ↔ I:ABU-N"),
    # 15-16 S156 ↔ ABU-N
    ("chainID A and resid 156 and name O",   "chainID H and resname ABU and name N",
     "A:Ser156-O ↔ H:ABU-N"),
    ("chainID C and resid 156 and name O",   "chainID I and resname ABU and name N",
     "C:Ser156-O ↔ I:ABU-N"),

    # R207 trio
    # 17-18 E153-R207
    ("chainID A and resid 153 and name CD", "chainID A and resid 207 and name NH1",
     "A:Glu153-CD ↔ A:Arg207-NH1"),
    ("chainID C and resid 153 and name CD", "chainID C and resid 207 and name NH1",
     "C:Glu153-CD ↔ C:Arg207-NH1"),
    # 19-20 E155-R207
    ("chainID A and resid 155 and name CD", "chainID A and resid 207 and name NH2",
     "A:Glu155-CD ↔ A:Arg207-NH2"),
    ("chainID C and resid 155 and name CD", "chainID C and resid 207 and name NH2",
     "C:Glu155-CD ↔ C:Arg207-NH2"),

    # B9-B10 strands h-bonds
    # 21-24 207-196
    ("chainID A and resid 207 and name N", "chainID A and resid 196 and name O",
     "A:207-N ↔ A:196-O"),
    ("chainID C and resid 207 and name N", "chainID C and resid 196 and name O",
     "C:207-N ↔ C:196-O"),
    ("chainID A and resid 207 and name O", "chainID A and resid 196 and name N",
     "A:207-O ↔ A:196-N"),
    ("chainID C and resid 207 and name O", "chainID C and resid 196 and name N",
     "C:207-O ↔ C:196-N"),
    # 25-28 205-198
    ("chainID A and resid 205 and name N", "chainID A and resid 198 and name O",
     "A:205-N ↔ A:198-O"),
    ("chainID C and resid 205 and name N", "chainID C and resid 198 and name O",
     "C:205-N ↔ C:198-O"),
    ("chainID A and resid 205 and name O", "chainID A and resid 198 and name N",
     "A:205-O ↔ A:198-N"),
    ("chainID C and resid 205 and name O", "chainID C and resid 198 and name N",
     "C:205-O ↔ C:198-N"),
    # 29-32 203-200
    ("chainID A and resid 203 and name N", "chainID A and resid 200 and name O",
     "A:203-N ↔ A:200-O"),
    ("chainID C and resid 203 and name N", "chainID C and resid 200 and name O",
     "C:203-N ↔ C:200-O"),
    ("chainID A and resid 203 and name O", "chainID A and resid 200 and name N",
     "A:203-O ↔ A:200-N"),
    ("chainID C and resid 203 and name O", "chainID C and resid 200 and name N",
     "C:203-O ↔ C:200-N"),

]

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="Atom-pair distances from GROMACS trajectory")
p.add_argument("-s", "--tpr", required=True, help="Topology file (.tpr/.gro/.pdb)")
p.add_argument("-f", "--xtc", required=True, help="Trajectory file (.xtc/.trr)")
p.add_argument("-b", "--begin", type=float, default=None, help="Start time in ns (default: first frame)")
p.add_argument("-e", "--end",   type=float, default=None, help="End time in ns (default: last frame)")
p.add_argument("-dt", "--step", type=int,   default=None, help="Process every Nth frame (default: all)")
p.add_argument("-w", "--window", type=int, default=50, help="Smoothing window in frames (default: 50, 0=off)")
p.add_argument("-o", "--prefix", default="distances", help="Output file prefix")
args = p.parse_args()

# ── Load universe ────────────────────────────────────────────────────
u = mda.Universe(args.tpr, args.xtc)

# Convert ns → ps for MDAnalysis slicing
start_ps = args.begin * 1000 if args.begin is not None else None
stop_ps  = args.end   * 1000 if args.end   is not None else None
traj_slice = u.trajectory[:]  # default: all frames
if start_ps is not None or stop_ps is not None or args.step is not None:
    # Build start/stop/step for frame slicing via time
    start_frame = 0
    stop_frame  = u.trajectory.n_frames
    step_frame  = args.step or 1
    for i, ts in enumerate(u.trajectory):
        if start_ps is not None and ts.time < start_ps:
            start_frame = i + 1
        if stop_ps is not None and ts.time > stop_ps:
            stop_frame = i
            break
    traj_slice = u.trajectory[start_frame:stop_frame:step_frame]

print(f"Loaded: {u.trajectory.n_frames} total frames, {u.atoms.n_atoms} atoms")
print(f"Analysing: {len(traj_slice)} frames"
      f" ({traj_slice[0].time/1000:.2f}–{traj_slice[-1].time/1000:.2f} ns, step={args.step or 1})")

# ── Validate selections ─────────────────────────────────────────────
valid_pairs = []
for sel1, sel2, label in PAIRS:
    ag1, ag2 = u.select_atoms(sel1), u.select_atoms(sel2)
    if ag1.n_atoms == 0 or ag2.n_atoms == 0:
        print(f"  SKIP  {label}  (sel1: {ag1.n_atoms}, sel2: {ag2.n_atoms} atoms)")
        continue
    if ag1.n_atoms > 1 or ag2.n_atoms > 1:
        print(f"  WARN  {label}  using first atom (sel1: {ag1.n_atoms}, sel2: {ag2.n_atoms})")
    valid_pairs.append((sel1, sel2, label))
    print(f"  OK    {label}")

# ── Calculate distances ──────────────────────────────────────────────
# Pre-select atoms ONCE (avoid re-parsing selections every frame)
atom_pairs = []
for sel1, sel2, label in valid_pairs:
    a1 = u.select_atoms(sel1)[0]  # single Atom object
    a2 = u.select_atoms(sel2)[0]
    atom_pairs.append((a1.index, a2.index))

idx1 = np.array([p[0] for p in atom_pairs])
idx2 = np.array([p[1] for p in atom_pairs])

n_frames = len(traj_slice)
times = np.zeros(n_frames)
all_dists = np.zeros((len(valid_pairs), n_frames))

# Single loop, vectorized distance for ALL pairs at once
for i, ts in enumerate(traj_slice):
    times[i] = ts.time / 1000.0  # ps → ns
    pos = ts.positions  # direct reference, no copy
    diff = pos[idx1] - pos[idx2]  # (n_pairs, 3) in one shot
    all_dists[:, i] = np.sqrt((diff * diff).sum(axis=1))
    if (i + 1) % 2000 == 0 or i == n_frames - 1:
        print(f"  Frame {i+1}/{n_frames}")

# ── Save CSV ─────────────────────────────────────────────────────────
labels = [lbl for _, _, lbl in valid_pairs]
header = "time_ns," + ",".join(f'"{l}"' for l in labels)
out_data = np.column_stack([times] + [all_dists[j] for j in range(len(valid_pairs))])
np.savetxt(f"{args.prefix}.csv", out_data, delimiter=",", header=header, comments="")
print(f"Saved {args.prefix}.csv")

# ── Interactive Plotly figure ────────────────────────────────────────
# Group pairs by similar interactions for subplots

def smooth(y, window=50):
    """Running average with edge handling."""
    kernel = np.ones(window) / window
    return np.convolve(y, kernel, mode="same")

fig = make_subplots(rows=4, cols=1, shared_xaxes=True, vertical_spacing=0.06,
                    subplot_titles=(
                        "loop C cap",
                        "GABA binding",
                        "R207 trio",
                        "loop C strands h-bonds"
                    ))

# Assign pairs to subplot rows (2 traces each)
row_map = [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,]
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

for j, (_, _, label) in enumerate(valid_pairs):
    row = row_map[j] if j < len(row_map) else 3
    ref_d = all_dists[j, 0]
    raw = all_dists[j] - ref_d
    color = colors[j % len(colors)]

    # Raw trace (thin, semi-transparent)
    fig.add_trace(
        go.Scatter(
            x=times, y=raw, mode="lines", name=label,
            line=dict(width=0.5, color=color), opacity=0.3,
            hovertemplate=f"{label}<br>t=%{{x:.2f}} ns<br>d=%{{y:.2f}} Å<extra></extra>",
            legendgroup=label,
        ),
        row=row, col=1,
    )
    # Smoothed trace (bold, on top)
    if args.window > 1:
        fig.add_trace(
            go.Scatter(
                x=times, y=smooth(raw, args.window), mode="lines",
                name=f"{label} (avg {args.window}f)", showlegend=False,
                line=dict(width=2, color=color),
                hovertemplate=f"{label} smoothed<br>t=%{{x:.2f}} ns<br>d=%{{y:.2f}} Å<extra></extra>",
                legendgroup=label,
            ),
            row=row, col=1,
        )
    # add initial-frame reference as dashed line (toggles with main trace)
    '''
    color = colors[j % len(colors)]
    fig.add_trace(
        go.Scatter(
            x=[times[0], times[-1]], y=[ref_d, ref_d], mode="lines+text",
            line=dict(width=1, dash="dot", color=color),
            opacity=0.5, showlegend=False, legendgroup=label,
            text=[f"t₀={ref_d:.2f} Å", ""], textposition="top left",
            textfont=dict(size=10, color=color),
            hoverinfo="skip",
        ),
        row=row, col=1,
    
    )
    '''

for r in range(1, 4):
    fig.update_yaxes(title_text="Distance (Å)", row=r, col=1)
fig.update_xaxes(title_text="Time (ns)", row=3, col=1)

fig.update_layout(
    height=900, template="plotly_white",
    legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="center", x=0.5,
                font=dict(size=9)),
    title="Atom-Pair Distances",
)

out_html = f"{args.prefix}_distance_plots.html"
fig.write_html(out_html, include_plotlyjs="cdn")
print(f"Saved {out_html}")
fig.show()