#!/usr/bin/env python3
"""Calculate atom-pair distances from one or more replicas of the same system
and save to a single CSV in long format (one row per replica × frame), with a
'replica' column added alongside 'time_ns'."""

import argparse
import csv
import numpy as np
import MDAnalysis as mda

# ── Atom pair definitions ────────────────────────────────────────────
# Each entry: (selection1, selection2, label, group)
PAIRS = [

    # BS1
    # glycan h-bond
    ('chainID A and resid 2 and name O', 'chainID A and resid 196 and name NZ',
     'A:glyc2-O ↔ A:Lys196-NZ', 'BS_1'),

    # trio
    # TODO: oxygen distances instead of CD
    ("chainID A and resid 153 and name CD", "chainID A and resid 207 and name NH1",
     "A:Glu153-CD ↔ A:Arg207-NH1", 'BS_1'),
    ("chainID A and resid 155 and name CD", "chainID A and resid 207 and name NH2",
     "A:Glu155-CD ↔ A:Arg207-NH2", 'BS_1'),

    # loop C vs F45
    ("chainID A and resid 205 and name CA", "chainID B and resid 46 and name CA",
     "A:X205-CA ↔ B:Phe46-CA", "BS_1"),
    ("chainID A and resid 204 and name CA", "chainID B and resid 46 and name CA",
     "A:X204-CA ↔ B:Phe46-CA", "BS_1"),
    ("chainID A and resid 203 and name CA", "chainID B and resid 46 and name CA",
     "A:X203-CA ↔ B:Phe46-CA", "BS_1"),
    ("chainID A and resid 202 and name CA", "chainID B and resid 46 and name CA",
     "A:X202-CA ↔ B:Phe46-CA", "BS_1"),
    ("chainID A and resid 201 and name CA", "chainID B and resid 46 and name CA",
     "A:X201-CA ↔ B:Phe46-CA", "BS_1"),
    ("chainID A and resid 200 and name CA", "chainID B and resid 46 and name CA",
     "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    ("chainID A and resid 199 and name CA", "chainID B and resid 46 and name CA",
     "A:X199-CA ↔ B:Phe46-CA", "BS_1"),
    ("chainID A and resid 198 and name CA", "chainID B and resid 46 and name CA",
     "A:X198-CA ↔ B:Phe46-CA", "BS_1"),
    ("chainID A and resid 197 and name CA", "chainID B and resid 46 and name CA",
     "A:X197-CA ↔ B:Phe46-CA", "BS_1"),
    ("chainID A and resid 196 and name CA", "chainID B and resid 46 and name CA",
     "A:Lys196-CA ↔ B:Phe46-CA", "BS_1"),
    ("chainID A and resid 195 and name CA", "chainID B and resid 46 and name CA",
     "A:X195-CA ↔ B:Phe46-CA", "BS_1"),
    ("chainID A and resid 194 and name CA", "chainID B and resid 46 and name CA",
     "A:X194-CA ↔ B:Phe46-CA", "BS_1"),

    # K196 and K196 s-bridges
    ("chainID A and resid 197 and name NZ", "chainID A and resid 165 and name OE1",
     "A:Lys197-NZ ↔ A:Glu165-OE1", "BS_1"),
    ("chainID A and resid 197 and name NZ", "chainID A and resid 165 and name OE2",
     "A:Lys197-NZ ↔ A:Glu165-OE2", "BS_1"),
    ("chainID A and resid 196 and name NZ", "chainID A and resid 153 and name OE1",
     "A:Lys196-NZ ↔ A:Glu153-OE1", "BS_1"),
    ("chainID A and resid 196 and name NZ", "chainID A and resid 153 and name OE2",
     "A:Lys196-NZ ↔ A:Glu153-OE2", "BS_1"),

    # interfaces
    #TODO:  add E153/D184 (loop F)
    ("chainID A and resid 205 and name CA", "chainID B and resid 67 and name CA",
     "A:Tyr205-CA ↔ B:Arg67-CA", "BS_1"),
    ("chainID A and resid 31 and name CA", "chainID B and resid 15 and name CA",
     "A:Phe31-CA ↔ B:Phe15-CA", "BS_1"),
    ("chainID A and resid 157 and name CA", "chainID B and resid 117 and name CA",
     "A:Tyr157-CA ↔ B:Lys117-CA", "BS_1"),
    ("chainID A and resid 136 and name CA", "chainID B and resid 186 and name CA",
     "A:Cys136-CA ↔ B:Ser186-CA", "BS_1"),

    # h-bond loop C cap
    ("chainID A and resid 205 and name OH", "chainID A and resid 202 and name OG1",
     "A:Tyr205-OH ↔ A:Thr202-OG1", "BS_1"),
    ("chainID A and resid 162 and name CG", "chainID A and resid 202 and name OG1",
     "A:Asp162-CG ↔ A:Thr202-OG1", "BS_1"),


    ("chainID A and resid 205 and name OH", "chainID A and resid 156 and name O",
     "A:Tyr205-OH ↔ A:Ser156-O", "BS_1"),
    ("chainID A and resid 205 and name OH", "chainID A and resid 159 and name O",
     "A:Tyr205-OH ↔ A:Tyr159-O", "BS_1"),
    ("chainID A and resid 205 and name OH", "chainID B and resid 120 and name CZ",
     "A:Tyr205-OH ↔ B:Arg120-CZ", "BS_1"),



    # GABA contacts
    ("chainID A and resid 202 and name OG1", "chainID H and resname ABU and name C4",
     "A:Thr202-OG1 ↔ H:ABU-C4", "BS_1"),
    ("chainID B and resid 67 and name CZ", "chainID H and resname ABU and name C4",
     "B:Arg67-CZ ↔ H:ABU-C4", "BS_1"),
    ("chainID B and resid 130 and name OG1", "chainID H and resname ABU and name C4",
     "B:Thr130-OG1 ↔ H:ABU-C4", "BS_1"),
    ("chainID A and resid 155 and name CD", "chainID H and resname ABU and name N",
     "A:Glu155-CD ↔ H:ABU-N", "BS_1"),
    ("chainID A and resid 97 and name OH", "chainID H and resname ABU and name N",
     "A:Tyr97-OH ↔ H:ABU-N", "BS_1"),
    ("chainID A and resid 156 and name O", "chainID H and resname ABU and name N",
     "A:Ser156-O ↔ H:ABU-N", "BS_1"),


    # BS2

    ('chainID C and resid 2 and name O', 'chainID C and resid 196 and name NZ',
     'C:glyc2-O ↔ C:Lys197-NZ', 'BS_2'),

    # trio
    # TODO: oxygen distances instead of CD
    ("chainID C and resid 153 and name CD", "chainID C and resid 207 and name NH1",
     "C:Glu153-CD ↔ C:Arg207-NH1", 'BS_2'),
    ("chainID C and resid 155 and name CD", "chainID C and resid 207 and name NH2",
     "C:Glu155-CD ↔ C:Arg207-NH2", 'BS_2'),

    # loop C vs F45
    ("chainID C and resid 205 and name CA", "chainID D and resid 46 and name CA",
     "C:X205-CA ↔ D:Phe46-CA", "BS_2"),
    ("chainID C and resid 204 and name CA", "chainID D and resid 46 and name CA",
     "C:X204-CA ↔ D:Phe46-CA", "BS_2"),
    ("chainID C and resid 203 and name CA", "chainID D and resid 46 and name CA",
     "C:X203-CA ↔ D:Phe46-CA", "BS_2"),
    ("chainID C and resid 202 and name CA", "chainID D and resid 46 and name CA",
     "C:X202-CA ↔ D:Phe46-CA", "BS_2"),
    ("chainID C and resid 201 and name CA", "chainID D and resid 46 and name CA",
     "C:X201-CA ↔ D:Phe46-CA", "BS_2"),
    ("chainID C and resid 200 and name CA", "chainID D and resid 46 and name CA",
     "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),
    ("chainID C and resid 199 and name CA", "chainID D and resid 46 and name CA",
     "C:X199-CA ↔ D:Phe46-CA", "BS_2"),
    ("chainID C and resid 198 and name CA", "chainID D and resid 46 and name CA",
     "C:X198-CA ↔ D:Phe46-CA", "BS_2"),
    ("chainID C and resid 197 and name CA", "chainID D and resid 46 and name CA",
     "C:X197-CA ↔ D:Phe46-CA", "BS_2"),
    ("chainID C and resid 196 and name CA", "chainID D and resid 46 and name CA",
     "C:Lys196-CA ↔ D:Phe46-CA", "BS_2"),
    ("chainID C and resid 195 and name CA", "chainID D and resid 46 and name CA",
     "C:X195-CA ↔ D:Phe46-CA", "BS_2"),
    ("chainID C and resid 194 and name CA", "chainID D and resid 46 and name CA",
     "C:X194-CA ↔ D:Phe46-CA", "BS_2"),

    # K196 and K196 s-bridges
    ("chainID C and resid 197 and name NZ", "chainID C and resid 165 and name OE1",
     "C:Lys197-NZ ↔ C:Glu165-OE1", "BS_2"),
    ("chainID C and resid 197 and name NZ", "chainID C and resid 165 and name OE2",
     "C:Lys197-NZ ↔ C:Glu165-OE2", "BS_2"),
    ("chainID C and resid 196 and name NZ", "chainID C and resid 153 and name OE1",
     "C:Lys196-NZ ↔ C:Glu153-OE1", "BS_2"),
    ("chainID C and resid 196 and name NZ", "chainID C and resid 153 and name OE2",
     "C:Lys196-NZ ↔ C:Glu153-OE2", "BS_2"),

    # interfaces
    # TODO:  add E153/D184 (loop F)
    ("chainID C and resid 205 and name CA", "chainID D and resid 67 and name CA",
     "C:Tyr205-CA ↔ D:Arg67-CA", "BS_2"),
    ("chainID C and resid 31 and name CA", "chainID D and resid 15 and name CA",
     "C:Phe31-CA ↔ D:Phe15-CA", "BS_2"),
    ("chainID C and resid 157 and name CA", "chainID D and resid 117 and name CA",
     "C:Tyr157-CA ↔ D:Lys117-CA", "BS_2"),
    ("chainID C and resid 136 and name CA", "chainID D and resid 186 and name CA",
     "C:Cys136-CA ↔ D:Ser186-CA", "BS_2"),

    # h-bond loop C cap
    ("chainID C and resid 205 and name OH", "chainID C and resid 202 and name OG1",
     "C:Tyr205-OH ↔ C:Thr202-OG1", "BS_2"),
    ("chainID C and resid 162 and name CG", "chainID C and resid 202 and name OG1",
     "C:Asp162-CG ↔ C:Thr202-OG1", "BS_2"),

    ("chainID C and resid 205 and name OH", "chainID C and resid 156 and name O",
     "C:Tyr205-OH ↔ C:Ser156-O", "BS_2"),
    ("chainID C and resid 205 and name OH", "chainID C and resid 159 and name O",
     "C:Tyr205-OH ↔ C:Tyr159-O", "BS_2"),
    ("chainID C and resid 205 and name OH", "chainID D and resid 120 and name CZ",
     "C:Tyr205-OH ↔ D:Arg120-CZ", "BS_2"),

    # GABA contacts
    ("chainID C and resid 202 and name OG1", "chainID I and resname ABU and name C4",
     "C:Thr202-OG1 ↔ I:ABU-C4", "BS_2"),
    ("chainID D and resid 67 and name CZ", "chainID I and resname ABU and name C4",
     "D:Arg67-CZ ↔ I:ABU-C4", "BS_2"),
    ("chainID D and resid 130 and name OG1", "chainID I and resname ABU and name C4",
     "D:Thr130-OG1 ↔ I:ABU-C4", "BS_2"),
    ("chainID C and resid 155 and name CD", "chainID I and resname ABU and name N",
     "C:Glu155-CD ↔ I:ABU-N", "BS_2"),
    ("chainID C and resid 97 and name OH", "chainID I and resname ABU and name N",
     "C:Tyr97-OH ↔ I:ABU-N", "BS_2"),
    ("chainID C and resid 156 and name O", "chainID I and resname ABU and name N",
     "C:Ser156-O ↔ I:ABU-N", "BS_2"),

    # interface a/b
    ("chainID B and resid 205 and name CA", "chainID C and resid 43 and name CA",
     "B:Ser205-CA ↔ C:Asp43-CA", "a/b"),

    # interface a/g
    ("chainID D and resid 205 and name CA", "chainID E and resid 58 and name CA",
     "D:Ser205-CA ↔ E:Tyr58-CA", "a/g"),

    # interface g/a
    ("chainID E and resid 215 and name CA", "chainID A and resid 43 and name CA",
     "E:Phe215-CA ↔ A:Asp43-CA", "g/a"),

]

PAIRS_UNUSED = [
    # BS_opening multi-sub
    ("chainID A and resid 200 and name CA", "chainID B and resid 46 and name CA",
     "A:Phe200-CA ↔ B:Phe46-CA", "BS_opening"),
    ("chainID B and resid 205 and name CA", "chainID C and resid 43 and name CA",
     "B:Ser205-CA ↔ C:Asp43-CA", "BS_opening"),
    ("chainID C and resid 200 and name CA", "chainID D and resid 46 and name CA",
     "C:Phe200-CA ↔ D:Phe46-CA", "BS_opening"),
    ("chainID D and resid 205 and name CA", "chainID E and resid 58 and name CA",
     "D:Ser205-CA ↔ E:Tyr58-CA", "BS_opening"),
    ("chainID E and resid 215 and name CA", "chainID A and resid 43 and name CA",
     "E:Phe215-CA ↔ A:Asp43-CA", "BS_opening"),

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
p = argparse.ArgumentParser(description="Atom-pair distances from GROMACS trajectories (multi-replica)")
p.add_argument("-s", "--tpr", required=True, help="Topology file (.tpr/.gro/.pdb), shared across replicas")
p.add_argument("-f", "--xtc", required=True, nargs="+",
               help="Trajectory file(s); one per replica (labelled replica1, replica2, …)")
p.add_argument("-b", "--begin", type=float, default=None, nargs="+",
               help="Start time(s) in ns: one value for all replicas, or one per replica")
p.add_argument("-e", "--end", type=float, default=None, nargs="+",
               help="End time(s) in ns: one value for all replicas, or one per replica")
p.add_argument("-dt", "--step", type=int, default=None, nargs="+",
               help="Frame step(s): one value for all replicas, or one per replica")
p.add_argument("-o", "--prefix", default="distances", help="Output file prefix")
args = p.parse_args()


def per_replica(arg, n):
    """Expand a 1- or n-list to length n; None → list of Nones."""
    if arg is None:
        return [None] * n
    if len(arg) == 1:
        return list(arg) * n
    if len(arg) == n:
        return list(arg)
    raise SystemExit(f"Expected 1 or {n} values, got {len(arg)}")


n_replicas = len(args.xtc)
begins = per_replica(args.begin, n_replicas)
ends   = per_replica(args.end,   n_replicas)
steps  = per_replica(args.step,  n_replicas)

# ── Validate selections (once, against first replica's topology) ────
u0 = mda.Universe(args.tpr, args.xtc[0])
print(f"Topology: {u0.atoms.n_atoms} atoms; {n_replicas} replica(s)")

valid_pairs = []
for sel1, sel2, label, group in PAIRS:
    ag1, ag2 = u0.select_atoms(sel1), u0.select_atoms(sel2)
    if ag1.n_atoms == 0 or ag2.n_atoms == 0:
        print(f"  SKIP  {label}  (sel1: {ag1.n_atoms}, sel2: {ag2.n_atoms} atoms)")
        continue
    if ag1.n_atoms > 1 or ag2.n_atoms > 1:
        print(f"  WARN  {label}  using first atom (sel1: {ag1.n_atoms}, sel2: {ag2.n_atoms})")
    valid_pairs.append((sel1, sel2, label, group))
    print(f"  OK    [{group}]  {label}")

# Pre-resolve atom indices once — topology is shared, so indices are stable
idx1 = np.array([u0.select_atoms(s1)[0].index for s1, _, _, _ in valid_pairs])
idx2 = np.array([u0.select_atoms(s2)[0].index for _, s2, _, _ in valid_pairs])

# ── Process each replica ────────────────────────────────────────────
replica_results = []  # list of (rep_label, times, all_dists)

for r_idx, xtc in enumerate(args.xtc):
    rep_label = f"replica{r_idx + 1}"
    print(f"\n── {rep_label}: {xtc} ──")

    u = u0 if r_idx == 0 else mda.Universe(args.tpr, xtc)

    begin, end, step = begins[r_idx], ends[r_idx], steps[r_idx]
    start_ps = begin * 1000 if begin is not None else None
    stop_ps  = end   * 1000 if end   is not None else None
    traj_slice = u.trajectory[:]
    if start_ps is not None or stop_ps is not None or step is not None:
        start_frame = 0
        stop_frame  = u.trajectory.n_frames
        step_frame  = step or 1
        for i, ts in enumerate(u.trajectory):
            if start_ps is not None and ts.time < start_ps:
                start_frame = i + 1
            if stop_ps is not None and ts.time > stop_ps:
                stop_frame = i
                break
        traj_slice = u.trajectory[start_frame:stop_frame:step_frame]

    print(f"  Total {u.trajectory.n_frames} frames; analysing {len(traj_slice)}"
          f" ({traj_slice[0].time / 1000:.2f}–{traj_slice[-1].time / 1000:.2f} ns, step={step or 1})")

    n_frames = len(traj_slice)
    times = np.zeros(n_frames)
    all_dists = np.zeros((len(valid_pairs), n_frames))

    for i, ts in enumerate(traj_slice):
        times[i] = ts.time / 1000.0  # ps → ns
        pos = ts.positions
        diff = pos[idx1] - pos[idx2]
        all_dists[:, i] = np.sqrt((diff * diff).sum(axis=1))
        if (i + 1) % 2000 == 0 or i == n_frames - 1:
            print(f"  Frame {i + 1}/{n_frames}")

    replica_results.append((rep_label, times, all_dists))

# ── Save CSV (long format) ──────────────────────────────────────────
# Columns: time_ns, replica, "group|label1", "group|label2", …
# csv.writer handles the mixed string/float types cleanly.
header = ["time_ns", "replica"] + [f"{g}|{l}" for _, _, l, g in valid_pairs]
out_path = f"{args.prefix}_distances.csv"
with open(out_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(header)
    for rep_label, times, all_dists in replica_results:
        for i in range(len(times)):
            row = [f"{times[i]:.4f}", rep_label] + [f"{all_dists[j, i]:.4f}" for j in range(all_dists.shape[0])]
            w.writerow(row)

total_rows = sum(len(t) for _, t, _ in replica_results)
print(f"\nSaved {out_path} ({total_rows} rows across {n_replicas} replicas)")