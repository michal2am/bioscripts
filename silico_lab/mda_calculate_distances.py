#!/usr/bin/env python3
"""Calculate distances between selected atom pairs from a GROMACS trajectory and save to CSV."""

import argparse
import numpy as np
import MDAnalysis as mda

# ── Atom pair definitions ────────────────────────────────────────────
# Each entry: (selection1, selection2, label, group)
PAIRS = [

    # loop C
    # 1-2 F200 ↔ T202 h-bond cap
    ("chainID A and resid 205 and name OH", "chainID A and resid 202 and name OG1",
     "A:Tyr205-OH ↔ A:Thr202-OG1", "loop C"),
    ("chainID C and resid 205 and name OH", "chainID C and resid 202 and name OG1",
     "C:Tyr205-OH ↔ C:Thr202-OG1", "loop C"),
    # 3-4 F200 ↔ F46 BS opening
    ("chainID A and resid 200 and name CG", "chainID B and resid 46 and name CG",
     "A:200-CG ↔ B:46-CG", "loop C"),
    ("chainID C and resid 200 and name CG", "chainID D and resid 46 and name CG",
     "C:200-CG ↔ D:46-CG", "loop C"),

    # GABA electro
    # 5-6 T202 ↔ ABU-C4
    ("chainID A and resid 202 and name OG1", "chainID H and resname ABU and name C4",
     "A:Thr202-OG1 ↔ H:ABU-C4", "GABA electro"),
    ("chainID C and resid 202 and name OG1", "chainID I and resname ABU and name C4",
     "C:Thr202-OG1 ↔ I:ABU-C4", "GABA electro"),
    # 7-8 R67 ↔ ABU-C4
    ("chainID B and resid 67 and name CZ", "chainID H and resname ABU and name C4",
     "B:Arg67-CZ ↔ H:ABU-C4", "GABA electro"),
    ("chainID D and resid 67 and name CZ", "chainID I and resname ABU and name C4",
     "D:Arg67-CZ ↔ I:ABU-C4", "GABA electro"),
    # 9-10 T130 ↔ ABU-C4
    ("chainID B and resid 130 and name OG1", "chainID H and resname ABU and name C4",
     "B:Thr130-OG1 ↔ H:ABU-C4", "GABA electro"),
    ("chainID D and resid 130 and name OG1", "chainID I and resname ABU and name C4",
     "D:Thr130-OG1 ↔ I:ABU-C4", "GABA electro"),
    # 11-12 E155 ↔ ABU-N
    ("chainID A and resid 155 and name CD", "chainID H and resname ABU and name N",
     "A:Glu155-CD ↔ H:ABU-N", "GABA electro"),
    ("chainID C and resid 155 and name CD", "chainID I and resname ABU and name N",
     "C:Glu155-CD ↔ I:ABU-N", "GABA electro"),
    # 13-14 Y97 ↔ ABU-N
    ("chainID A and resid 97 and name OH", "chainID H and resname ABU and name N",
     "A:Tyr97-OH ↔ H:ABU-N", "GABA electro"),
    ("chainID C and resid 97 and name OH", "chainID I and resname ABU and name N",
     "C:Tyr97-OH ↔ I:ABU-N", "GABA electro"),
    # 15-16 S156 ↔ ABU-N
    ("chainID A and resid 156 and name O", "chainID H and resname ABU and name N",
     "A:Ser156-O ↔ H:ABU-N", "GABA electro"),
    ("chainID C and resid 156 and name O", "chainID I and resname ABU and name N",
     "C:Ser156-O ↔ I:ABU-N", "GABA electro"),

    # R207 trio
    # 17-18 E153-R207
    ("chainID A and resid 153 and name CD", "chainID A and resid 207 and name NH1",
     "A:Glu153-CD ↔ A:Arg207-NH1", "R207 trio"),
    ("chainID C and resid 153 and name CD", "chainID C and resid 207 and name NH1",
     "C:Glu153-CD ↔ C:Arg207-NH1", "R207 trio"),
    # 19-20 E155-R207
    ("chainID A and resid 155 and name CD", "chainID A and resid 207 and name NH2",
     "A:Glu155-CD ↔ A:Arg207-NH2", "R207 trio"),
    ("chainID C and resid 155 and name CD", "chainID C and resid 207 and name NH2",
     "C:Glu155-CD ↔ C:Arg207-NH2", "R207 trio"),

    # B9-B10 strands h-bonds
    # 21-24 207-196
    ("chainID A and resid 207 and name N", "chainID A and resid 196 and name O",
     "A:207-N ↔ A:196-O", "B9-B10 strands h-bonds"),
    ("chainID C and resid 207 and name N", "chainID C and resid 196 and name O",
     "C:207-N ↔ C:196-O", "B9-B10 strands h-bonds"),
    ("chainID A and resid 207 and name O", "chainID A and resid 196 and name N",
     "A:207-O ↔ A:196-N", "B9-B10 strands h-bonds"),
    ("chainID C and resid 207 and name O", "chainID C and resid 196 and name N",
     "C:207-O ↔ C:196-N", "B9-B10 strands h-bonds"),
    # 25-28 205-198
    ("chainID A and resid 205 and name N", "chainID A and resid 198 and name O",
     "A:205-N ↔ A:198-O", "B9-B10 strands h-bonds"),
    ("chainID C and resid 205 and name N", "chainID C and resid 198 and name O",
     "C:205-N ↔ C:198-O", "B9-B10 strands h-bonds"),
    ("chainID A and resid 205 and name O", "chainID A and resid 198 and name N",
     "A:205-O ↔ A:198-N", "B9-B10 strands h-bonds"),
    ("chainID C and resid 205 and name O", "chainID C and resid 198 and name N",
     "C:205-O ↔ C:198-N", "B9-B10 strands h-bonds"),
    # 29-32 203-200
    ("chainID A and resid 203 and name N", "chainID A and resid 200 and name O",
     "A:203-N ↔ A:200-O", "B9-B10 strands h-bonds"),
    ("chainID C and resid 203 and name N", "chainID C and resid 200 and name O",
     "C:203-N ↔ C:200-O", "B9-B10 strands h-bonds"),
    ("chainID A and resid 203 and name O", "chainID A and resid 200 and name N",
     "A:203-O ↔ A:200-N", "B9-B10 strands h-bonds"),
    ("chainID C and resid 203 and name O", "chainID C and resid 200 and name N",
     "C:203-O ↔ C:200-N", "B9-B10 strands h-bonds"),

]

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="Atom-pair distances from GROMACS trajectory")
p.add_argument("-s", "--tpr", required=True, help="Topology file (.tpr/.gro/.pdb)")
p.add_argument("-f", "--xtc", required=True, help="Trajectory file (.xtc/.trr)")
p.add_argument("-b", "--begin", type=float, default=None, help="Start time in ns (default: first frame)")
p.add_argument("-e", "--end", type=float, default=None, help="End time in ns (default: last frame)")
p.add_argument("-dt", "--step", type=int, default=None, help="Process every Nth frame (default: all)")
p.add_argument("-o", "--prefix", default="distances", help="Output file prefix")
args = p.parse_args()

# ── Load universe ────────────────────────────────────────────────────
u = mda.Universe(args.tpr, args.xtc)

# Convert ns → ps for MDAnalysis slicing
start_ps = args.begin * 1000 if args.begin is not None else None
stop_ps = args.end * 1000 if args.end is not None else None
traj_slice = u.trajectory[:]  # default: all frames
if start_ps is not None or stop_ps is not None or args.step is not None:
    # Build start/stop/step for frame slicing via time
    start_frame = 0
    stop_frame = u.trajectory.n_frames
    step_frame = args.step or 1
    for i, ts in enumerate(u.trajectory):
        if start_ps is not None and ts.time < start_ps:
            start_frame = i + 1
        if stop_ps is not None and ts.time > stop_ps:
            stop_frame = i
            break
    traj_slice = u.trajectory[start_frame:stop_frame:step_frame]

print(f"Loaded: {u.trajectory.n_frames} total frames, {u.atoms.n_atoms} atoms")
print(f"Analysing: {len(traj_slice)} frames"
      f" ({traj_slice[0].time / 1000:.2f}–{traj_slice[-1].time / 1000:.2f} ns, step={args.step or 1})")

# ── Validate selections ─────────────────────────────────────────────
valid_pairs = []
for sel1, sel2, label, group in PAIRS:
    ag1, ag2 = u.select_atoms(sel1), u.select_atoms(sel2)
    if ag1.n_atoms == 0 or ag2.n_atoms == 0:
        print(f"  SKIP  {label}  (sel1: {ag1.n_atoms}, sel2: {ag2.n_atoms} atoms)")
        continue
    if ag1.n_atoms > 1 or ag2.n_atoms > 1:
        print(f"  WARN  {label}  using first atom (sel1: {ag1.n_atoms}, sel2: {ag2.n_atoms})")
    valid_pairs.append((sel1, sel2, label, group))
    print(f"  OK    [{group}]  {label}")

# ── Calculate distances ──────────────────────────────────────────────
# Pre-select atoms ONCE (avoid re-parsing selections every frame)
atom_pairs = []
for sel1, sel2, label, group in valid_pairs:
    a1 = u.select_atoms(sel1)[0]  # single Atom object
    a2 = u.select_atoms(sel2)[0]
    atom_pairs.append((a1.index, a2.index))

idx1 = np.array([pr[0] for pr in atom_pairs])
idx2 = np.array([pr[1] for pr in atom_pairs])

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
        print(f"  Frame {i + 1}/{n_frames}")

# ── Save CSV ─────────────────────────────────────────────────────────
# Encode group in column header as "group|label" so the plot script can
# assign each data series to a subplot panel by its group.
header = "time_ns," + ",".join(f'"{g}|{l}"' for _, _, l, g in valid_pairs)
out_data = np.column_stack([times] + [all_dists[j] for j in range(len(valid_pairs))])
np.savetxt(f"{args.prefix}_distances.csv", out_data, delimiter=",", header=header, comments="")
print(f"Saved {args.prefix}_distances.csv")