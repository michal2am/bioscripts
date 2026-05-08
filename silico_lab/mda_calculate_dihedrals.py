#!/usr/bin/env python3
"""Calculate backbone (phi, psi) and side-chain (chi1-chi4) dihedral angles
for selected residues from a GROMACS trajectory and save to CSV."""

import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals

# ── Residue definitions ──────────────────────────────────────────────
# Each entry: (selection, label) — selection should match a single residue
RESIDUES = [
    ("chainID A and resid 200", "A:200"),
    ("chainID C and resid 200", "C:200"),
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
p = argparse.ArgumentParser(description="Per-residue dihedral angles from GROMACS trajectory")
p.add_argument("-s", "--tpr", required=True, help="Topology file (.tpr/.gro/.pdb)")
p.add_argument("-f", "--xtc", required=True, help="Trajectory file (.xtc/.trr)")
p.add_argument("-b", "--begin", type=float, default=None, help="Start time in ns (default: first frame)")
p.add_argument("-e", "--end", type=float, default=None, help="End time in ns (default: last frame)")
p.add_argument("-dt", "--step", type=int, default=None, help="Process every Nth frame (default: all)")
p.add_argument("-o", "--prefix", default="dihedrals", help="Output file prefix")
args = p.parse_args()

# ── Load universe ────────────────────────────────────────────────────
u = mda.Universe(args.tpr, args.xtc)

# Convert ns → ps for MDAnalysis slicing
start_ps = args.begin * 1000 if args.begin is not None else None
stop_ps = args.end * 1000 if args.end is not None else None
traj_slice = u.trajectory[:]  # default: all frames
if start_ps is not None or stop_ps is not None or args.step is not None:
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

# ── Build dihedral list ─────────────────────────────────────────────
# Each entry: (label, group, [idx0, idx1, idx2, idx3])
dihedrals = []

for sel, res_label in RESIDUES:
    res_ag = u.select_atoms(sel)
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
        a = u.select_atoms(f"chainID {ch} and resid {rid} and name {name}")
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

# ── Calculate dihedrals ─────────────────────────────────────────────
idx0 = np.array([d[2][0] for d in dihedrals])
idx1 = np.array([d[2][1] for d in dihedrals])
idx2 = np.array([d[2][2] for d in dihedrals])
idx3 = np.array([d[2][3] for d in dihedrals])

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

# ── Save CSV ─────────────────────────────────────────────────────────
# Encode group in column header as "group|label" so the plot script can
# assign each data series to a subplot panel by its group.
header = "time_ns," + ",".join(f'"{g}|{l}"' for l, g, _ in dihedrals)
out_data = np.column_stack([times] + [all_angles[j] for j in range(len(dihedrals))])
np.savetxt(f"{args.prefix}_dihedrals.csv", out_data, delimiter=",", header=header, comments="")
print(f"Saved {args.prefix}_dihedrals.csv")