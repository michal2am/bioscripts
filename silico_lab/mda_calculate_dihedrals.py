#!/usr/bin/env python3
"""Calculate backbone (phi, psi) and side-chain (chi1-chi4) dihedral angles
for selected residues across one or more replicas of the same system, and
save to a single CSV in long format with a 'replica' column."""

import argparse
import csv
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals

# ── Residue definitions ──────────────────────────────────────────────
# Each entry: (selection, label) — selection should match a single residue
RESIDUES = [
    ("chainID A and resid 196", "A:196"),
    ("chainID A and resid 197", "A:197"),
    ("chainID C and resid 196", "C:196"),
    ("chainID C and resid 197", "C:197"),
]

# Standard chi-angle atom names per residue type
# (Gly and Ala have no chi angles)
CHI_ATOMS = {
    "ARG": [("N","CA","CB","CG"),  ("CA","CB","CG","CD"), ("CB","CG","CD","NE"), ("CG","CD","NE","CZ")],
    "ASN": [("N","CA","CB","CG"),  ("CA","CB","CG","OD1")],
    "ASP": [("N","CA","CB","CG"),  ("CA","CB","CG","OD1")],
    "CYS": [("N","CA","CB","SG")],
    "GLN": [("N","CA","CB","CG"),  ("CA","CB","CG","CD"), ("CB","CG","CD","OE1")],
    "GLU": [("N","CA","CB","CG"),  ("CA","CB","CG","CD"), ("CB","CG","CD","OE1")],
    "HIS": [("N","CA","CB","CG"),  ("CA","CB","CG","ND1")],
    "ILE": [("N","CA","CB","CG1"), ("CA","CB","CG1","CD1")],
    "LEU": [("N","CA","CB","CG"),  ("CA","CB","CG","CD1")],
    "LYS": [("N","CA","CB","CG"),  ("CA","CB","CG","CD"), ("CB","CG","CD","CE"), ("CG","CD","CE","NZ")],
    "MET": [("N","CA","CB","CG"),  ("CA","CB","CG","SD"), ("CB","CG","SD","CE")],
    "PHE": [("N","CA","CB","CG"),  ("CA","CB","CG","CD1")],
    "PRO": [("N","CA","CB","CG"),  ("CA","CB","CG","CD")],
    "SER": [("N","CA","CB","OG")],
    "THR": [("N","CA","CB","OG1")],
    "TRP": [("N","CA","CB","CG"),  ("CA","CB","CG","CD1")],
    "TYR": [("N","CA","CB","CG"),  ("CA","CB","CG","CD1")],
    "VAL": [("N","CA","CB","CG1")],
}

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="Per-residue dihedral angles from GROMACS trajectories (multi-replica)")
p.add_argument("-s", "--tpr", required=True, help="Topology file (.tpr/.gro/.pdb), shared across replicas")
p.add_argument("-f", "--xtc", required=True, nargs="+",
               help="Trajectory file(s); one per replica (labelled replica1, replica2, …)")
p.add_argument("-b", "--begin", type=float, default=None, nargs="+",
               help="Start time(s) in ns: one value for all replicas, or one per replica")
p.add_argument("-e", "--end", type=float, default=None, nargs="+",
               help="End time(s) in ns: one value for all replicas, or one per replica")
p.add_argument("-dt", "--step", type=int, default=None, nargs="+",
               help="Frame step(s): one value for all replicas, or one per replica")
p.add_argument("-o", "--prefix", default="dihedrals", help="Output file prefix")
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

# ── Resolve dihedrals once on first replica's topology ──────────────
# Each entry: (label, group, [idx0, idx1, idx2, idx3])
u0 = mda.Universe(args.tpr, args.xtc[0])
print(f"Topology: {u0.atoms.n_atoms} atoms; {n_replicas} replica(s)")

dihedrals = []
for sel, res_label in RESIDUES:
    res_ag = u0.select_atoms(sel)
    if res_ag.n_atoms == 0:
        print(f"  SKIP  {res_label}  (empty selection)")
        continue
    res = res_ag.residues[0]
    resname = res.resname
    resid = res.resid
    chain = res.atoms[0].chainID  # assumes chainID is set (matches selection syntax)
    group = f"{res_label} {resname}"
    print(f"  Residue [{group}] resid={resid} chain={chain}")

    def aidx(name, rid=resid, ch=chain):
        a = u0.select_atoms(f"chainID {ch} and resid {rid} and name {name}")
        return a[0].index if a.n_atoms == 1 else None

    # phi: prev.C – N – CA – C
    indices = [aidx("C", resid - 1), aidx("N"), aidx("CA"), aidx("C")]
    if all(i is not None for i in indices):
        dihedrals.append((f"{res_label}:phi", group, indices))
        print(f"    OK    phi")
    else:
        print(f"    SKIP  phi  (missing atoms)")

    # psi: N – CA – C – next.N
    indices = [aidx("N"), aidx("CA"), aidx("C"), aidx("N", resid + 1)]
    if all(i is not None for i in indices):
        dihedrals.append((f"{res_label}:psi", group, indices))
        print(f"    OK    psi")
    else:
        print(f"    SKIP  psi  (missing atoms)")

    # chi1, chi2, … by residue type
    if resname in CHI_ATOMS:
        for k, atom_names in enumerate(CHI_ATOMS[resname], start=1):
            indices = [aidx(n) for n in atom_names]
            if all(i is not None for i in indices):
                dihedrals.append((f"{res_label}:chi{k}", group, indices))
                print(f"    OK    chi{k}  ({'-'.join(atom_names)})")
            else:
                print(f"    SKIP  chi{k}  (missing atoms)")
    else:
        print(f"    note: no chi angles defined for {resname}")

# Stack indices into arrays once — topology is shared, so indices are stable
idx0 = np.array([d[2][0] for d in dihedrals])
idx1 = np.array([d[2][1] for d in dihedrals])
idx2 = np.array([d[2][2] for d in dihedrals])
idx3 = np.array([d[2][3] for d in dihedrals])

# ── Process each replica ────────────────────────────────────────────
replica_results = []  # list of (rep_label, times, all_angles)

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
    all_angles = np.zeros((len(dihedrals), n_frames))

    # Single loop, vectorized dihedral calculation for ALL angles at once
    for i, ts in enumerate(traj_slice):
        times[i] = ts.time / 1000.0  # ps → ns
        pos = ts.positions
        # calc_dihedrals returns radians → convert to degrees
        all_angles[:, i] = np.degrees(calc_dihedrals(pos[idx0], pos[idx1], pos[idx2], pos[idx3]))
        if (i + 1) % 2000 == 0 or i == n_frames - 1:
            print(f"  Frame {i + 1}/{n_frames}")

    replica_results.append((rep_label, times, all_angles))

# ── Save CSV (long format) ──────────────────────────────────────────
# Columns: time_ns, replica, "group|label1", "group|label2", …
# csv.writer handles the mixed string/float types cleanly.
header = ["time_ns", "replica"] + [f"{g}|{l}" for l, g, _ in dihedrals]
out_path = f"{args.prefix}_dihedrals.csv"
with open(out_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(header)
    for rep_label, times, all_angles in replica_results:
        for i in range(len(times)):
            row = [f"{times[i]:.4f}", rep_label] + [f"{all_angles[j, i]:.4f}" for j in range(all_angles.shape[0])]
            w.writerow(row)

total_rows = sum(len(t) for _, t, _ in replica_results)
print(f"\nSaved {out_path} ({total_rows} rows across {n_replicas} replicas)")