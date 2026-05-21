#!/usr/bin/env python3
"""
TICA on B9-B10 antiparallel beta strand H-bond distances + GABA contacts.
Joint TICA across one or more replicas; separate TICA per binding site.

Requires: numpy, pandas, MDAnalysis, deeptime
"""

import argparse
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_bonds, calc_dihedrals
from deeptime.decomposition import TICA


SITES = {
    "site1": [
        # GABA
        ("chainID A and resid 202 and name OG1", "chainID H and resname ABU and name C4", "Thr202-OG1↔ABU-C4"),
        ("chainID B and resid 67 and name CZ",   "chainID H and resname ABU and name C4", "Arg67-CZ↔ABU-C4"),
        ("chainID B and resid 130 and name OG1", "chainID H and resname ABU and name C4", "Thr130-OG1↔ABU-C4"),
        ("chainID A and resid 155 and name CD",  "chainID H and resname ABU and name N",  "Glu155-CD↔ABU-N"),
        ("chainID A and resid 97 and name OH",   "chainID H and resname ABU and name N",  "Tyr97-OH↔ABU-N"),
        ("chainID A and resid 156 and name O",   "chainID H and resname ABU and name N",  "Ser156-O↔ABU-N"),
        # h-bonds
        ("chainID A and resid 207 and name N", "chainID A and resid 196 and name O", "207-N↔196-O"),
        ("chainID A and resid 207 and name O", "chainID A and resid 196 and name N", "207-O↔196-N"),
        ("chainID A and resid 205 and name N", "chainID A and resid 198 and name O", "205-N↔198-O"),
        ("chainID A and resid 205 and name O", "chainID A and resid 198 and name N", "205-O↔198-N"),
        ("chainID A and resid 203 and name N", "chainID A and resid 200 and name O", "203-N↔200-O"),
        ("chainID A and resid 203 and name O", "chainID A and resid 200 and name N", "203-O↔200-N"),
        # trio
        ("chainID A and resid 153 and name CD", "chainID A and resid 207 and name NH1","153-CD↔207-NH1"),
        ("chainID A and resid 155 and name CD", "chainID A and resid 207 and name NH2","155-CD↔207-NH2"),
        # cap h-bond
        ("chainID A and resid 205 and name OH", "chainID A and resid 202 and name OG1","205-OH↔202-OG1"),
        # interfaces
        #("chainID A and resid 200 and name CG", "chainID B and resid 46 and name CG", "200-CG↔46-CG"),
    ],
    "site2": [
        # GABA
        ("chainID C and resid 202 and name OG1", "chainID I and resname ABU and name C4", "Thr202-OG1↔ABU-C4"),
        ("chainID D and resid 67 and name CZ",   "chainID I and resname ABU and name C4", "Arg67-CZ↔ABU-C4"),
        ("chainID D and resid 130 and name OG1", "chainID I and resname ABU and name C4", "Thr130-OG1↔ABU-C4"),
        ("chainID C and resid 155 and name CD",  "chainID I and resname ABU and name N",  "Glu155-CD↔ABU-N"),
        ("chainID C and resid 97 and name OH",   "chainID I and resname ABU and name N",  "Tyr97-OH↔ABU-N"),
        ("chainID C and resid 156 and name O",   "chainID I and resname ABU and name N",  "Ser156-O↔ABU-N"),
        # h-bonds
        ("chainID C and resid 207 and name N", "chainID C and resid 196 and name O", "207-N↔196-O"),
        ("chainID C and resid 207 and name O", "chainID C and resid 196 and name N", "207-O↔196-N"),
        ("chainID C and resid 205 and name N", "chainID C and resid 198 and name O", "205-N↔198-O"),
        ("chainID C and resid 205 and name O", "chainID C and resid 198 and name N", "205-O↔198-N"),
        ("chainID C and resid 203 and name N", "chainID C and resid 200 and name O", "203-N↔200-O"),
        ("chainID C and resid 203 and name O", "chainID C and resid 200 and name N", "203-O↔200-N"),
        # trio
        ("chainID C and resid 153 and name CD", "chainID C and resid 207 and name NH1","153-CD↔207-NH1"),
        ("chainID C and resid 155 and name CD", "chainID C and resid 207 and name NH2","155-CD↔207-NH2"),
        # cap h-bond
        ("chainID C and resid 205 and name OH", "chainID C and resid 202 and name OG1","205-OH↔202-OG1"),
        # interfaces
        #("chainID C and resid 200 and name CG", "chainID D and resid 46 and name CG", "200-CG↔46-CG"),
    ],
}

# Per-site auxiliary distances (computed separately, not used as TICA features)
AUX_DISTS = {
    "site1": ("chainID A and resid 200 and name CG", "chainID B and resid 46 and name CG", "200CG-46CG"),
    "site2": ("chainID C and resid 200 and name CG", "chainID D and resid 46 and name CG", "200CG-46CG"),
}

# Per-site chain containing the residues for dihedral computation
SITE_CHAIN = {"site1": "A", "site2": "C"}

# Backbone phi/psi residues
BACKBONE_RESIDS = [200, 201, 202, 203, 204, 205]

# Side chain chi residues (all available chi angles for each)
CHI_RESIDS = [200, 205]

# Standard chi angle atom definitions per amino acid (chi1, chi2, chi3, chi4)
CHI_ATOMS = {
    "ARG": [("N","CA","CB","CG"), ("CA","CB","CG","CD"), ("CB","CG","CD","NE"), ("CG","CD","NE","CZ")],
    "ASN": [("N","CA","CB","CG"), ("CA","CB","CG","OD1")],
    "ASP": [("N","CA","CB","CG"), ("CA","CB","CG","OD1")],
    "CYS": [("N","CA","CB","SG")],
    "GLN": [("N","CA","CB","CG"), ("CA","CB","CG","CD"), ("CB","CG","CD","OE1")],
    "GLU": [("N","CA","CB","CG"), ("CA","CB","CG","CD"), ("CB","CG","CD","OE1")],
    "HIS": [("N","CA","CB","CG"), ("CA","CB","CG","ND1")],
    "ILE": [("N","CA","CB","CG1"), ("CA","CB","CG1","CD1")],
    "LEU": [("N","CA","CB","CG"), ("CA","CB","CG","CD1")],
    "LYS": [("N","CA","CB","CG"), ("CA","CB","CG","CD"), ("CB","CG","CD","CE"), ("CG","CD","CE","NZ")],
    "MET": [("N","CA","CB","CG"), ("CA","CB","CG","SD"), ("CB","CG","SD","CE")],
    "PHE": [("N","CA","CB","CG"), ("CA","CB","CG","CD1")],
    "PRO": [("N","CA","CB","CG"), ("CA","CB","CG","CD")],
    "SER": [("N","CA","CB","OG")],
    "THR": [("N","CA","CB","OG1")],
    "TRP": [("N","CA","CB","CG"), ("CA","CB","CG","CD1")],
    "TYR": [("N","CA","CB","CG"), ("CA","CB","CG","CD1")],
    "VAL": [("N","CA","CB","CG1")],
}


p = argparse.ArgumentParser(description="Joint TICA across replicas, per binding site")
p.add_argument("-s", "--tpr", required=True, help="Topology file (.tpr/.gro/.pdb); shared across replicas")
p.add_argument("-f", "--xtc", required=True, nargs="+", help="Trajectory file(s); pass multiple for joint TICA")
p.add_argument("-b", "--begin", type=float, nargs="+", default=None,
               help="Start time(s) in ns: 1 value (all replicas) or 1 per replica (default: first frame)")
p.add_argument("-e", "--end", type=float, nargs="+", default=None,
               help="End time(s) in ns: 1 value (all replicas) or 1 per replica (default: last frame)")
p.add_argument("-dt", "--step", type=int, default=None, help="Process every Nth frame (default: all)")
p.add_argument("-lag", "--lag", type=float, default=1000.0, help="TICA lag time in ps (default: 1000)")
p.add_argument("-n", "--n-tics", type=int, default=None, help="Number of TICs to keep (default: all features)")
p.add_argument("-o", "--prefix", default="tica", help="Output file prefix")
args = p.parse_args()

step = args.step or 1

n_rep = len(args.xtc)

# Expand begin/end to per-replica lists (single value broadcasts; one-per-replica also accepted)
if args.begin is None:
    begins = [None] * n_rep
elif len(args.begin) == 1:
    begins = args.begin * n_rep
elif len(args.begin) == n_rep:
    begins = list(args.begin)
else:
    raise SystemExit(f"--begin: expected 1 or {n_rep} values, got {len(args.begin)}")

if args.end is None:
    ends = [None] * n_rep
elif len(args.end) == 1:
    ends = args.end * n_rep
elif len(args.end) == n_rep:
    ends = list(args.end)
else:
    raise SystemExit(f"--end: expected 1 or {n_rep} values, got {len(args.end)}")

# Resolve atom indices once from topology (shared across replicas)
u_ref = mda.Universe(args.tpr)
site_indices = {}
for site, pairs in SITES.items():
    idx1 = [u_ref.select_atoms(s1).indices[0] for s1, _, _ in pairs]
    idx2 = [u_ref.select_atoms(s2).indices[0] for _, s2, _ in pairs]
    site_indices[site] = (np.array(idx1), np.array(idx2))

# Resolve auxiliary distance atom indices (tracked alongside TICA but not used as features)
aux_indices = {}
for site, (s1, s2, label) in AUX_DISTS.items():
    aux_indices[site] = (
        u_ref.select_atoms(s1).indices[0],
        u_ref.select_atoms(s2).indices[0],
        label,
    )
aux_site_order = list(AUX_DISTS.keys())
aux_i1_arr = np.array([aux_indices[s][0] for s in aux_site_order])
aux_i2_arr = np.array([aux_indices[s][1] for s in aux_site_order])

# Resolve dihedral atom quads per site (backbone phi/psi + side chain chi)
site_dihed_quads = {}
site_dihed_labels = {}
for site, chain in SITE_CHAIN.items():
    quads = []
    labels = []
    # Backbone phi/psi
    for r in BACKBONE_RESIDS:
        # phi: C(i-1) - N(i) - CA(i) - C(i)
        quads.append([
            u_ref.select_atoms(f"chainID {chain} and resid {r-1} and name C").indices[0],
            u_ref.select_atoms(f"chainID {chain} and resid {r} and name N").indices[0],
            u_ref.select_atoms(f"chainID {chain} and resid {r} and name CA").indices[0],
            u_ref.select_atoms(f"chainID {chain} and resid {r} and name C").indices[0],
        ])
        labels.append(f"phi{r}")
        # psi: N(i) - CA(i) - C(i) - N(i+1)
        quads.append([
            u_ref.select_atoms(f"chainID {chain} and resid {r} and name N").indices[0],
            u_ref.select_atoms(f"chainID {chain} and resid {r} and name CA").indices[0],
            u_ref.select_atoms(f"chainID {chain} and resid {r} and name C").indices[0],
            u_ref.select_atoms(f"chainID {chain} and resid {r+1} and name N").indices[0],
        ])
        labels.append(f"psi{r}")
    # Side chain chi (all available per residue type)
    for r in CHI_RESIDS:
        resname = u_ref.select_atoms(f"chainID {chain} and resid {r}").residues[0].resname
        for ci, atoms in enumerate(CHI_ATOMS.get(resname, []), 1):
            quads.append([
                u_ref.select_atoms(f"chainID {chain} and resid {r} and name {a}").indices[0]
                for a in atoms
            ])
            labels.append(f"chi{ci}-{r}{resname}")
    site_dihed_quads[site] = np.array(quads)
    site_dihed_labels[site] = labels

n_dist = len(SITES["site1"])
n_dihed = len(site_dihed_labels["site1"])
print(f"Per-site features: {n_dist} distances + {n_dihed} dihedrals × 2 sin/cos = {n_dist + 2*n_dihed} total")
print(f"  site1 dihedrals: {', '.join(site_dihed_labels['site1'])}")
print(f"  site2 dihedrals: {', '.join(site_dihed_labels['site2'])}")
print(f"Aux distances (tracked, not in TICA): {', '.join(f'{s}|{aux_indices[s][2]}' for s in aux_site_order)}")
print()

# Collect features per replica (list of arrays per site for joint TICA fit)
features_per_rep = {site: [] for site in SITES}
aux_per_rep = {site: [] for site in AUX_DISTS}
times_per_rep = []
effective_dt_ps = None

for rep_idx, xtc in enumerate(args.xtc, 1):
    u = mda.Universe(args.tpr, xtc)
    if effective_dt_ps is None:
        effective_dt_ps = u.trajectory.dt * step
    begin = begins[rep_idx - 1]
    end = ends[rep_idx - 1]

    rep_times = []
    rep_feats = {site: [] for site in SITES}
    rep_aux = {site: [] for site in AUX_DISTS}
    for ts in u.trajectory:
        t_ns = ts.time / 1000.0
        if begin is not None and t_ns < begin:
            continue
        if end is not None and t_ns > end:
            break
        if (ts.frame % step) != 0:
            continue
        rep_times.append(t_ns)
        positions = u.atoms.positions
        for site, (i1, i2) in site_indices.items():
            d = calc_bonds(positions[i1], positions[i2], box=ts.dimensions)
            q = site_dihed_quads[site]
            angles = calc_dihedrals(
                positions[q[:, 0]], positions[q[:, 1]], positions[q[:, 2]], positions[q[:, 3]],
                box=ts.dimensions,
            )
            sc = np.empty(2 * len(angles))
            sc[0::2] = np.sin(angles)
            sc[1::2] = np.cos(angles)
            rep_feats[site].append(np.concatenate([d, sc]))
        d_aux = calc_bonds(positions[aux_i1_arr], positions[aux_i2_arr], box=ts.dimensions)
        for i, site in enumerate(aux_site_order):
            rep_aux[site].append(d_aux[i])

    for site in SITES:
        features_per_rep[site].append(np.array(rep_feats[site]))
    for site in AUX_DISTS:
        aux_per_rep[site].append(np.array(rep_aux[site]))
    times_per_rep.append(np.array(rep_times))
    print(f"Replica {rep_idx}: {len(rep_times)} frames from {xtc} (b={begin} ns, e={end} ns)")

lag_frames = max(1, int(round(args.lag / effective_dt_ps)))
print()
print(f"Lag time: {args.lag} ps = {lag_frames} frames (effective dt = {effective_dt_ps} ps)")
print(f"Total replicas: {len(args.xtc)}, total frames: {sum(len(t) for t in times_per_rep)}")
print()


# Fit TICA per site, joint across replicas
n_feat = n_dist + 2 * n_dihed
dim = args.n_tics if args.n_tics is not None else n_feat
dim = min(dim, n_feat)

projections = {}
eigenvalues = {}
timescales_ps = {}
loadings = {}

for site in SITES:
    tica = TICA(lagtime=lag_frames, dim=dim)
    # Pass list of arrays so lag-time covariance is accumulated within each replica only
    model = tica.fit(features_per_rep[site]).fetch_model()

    # Transform each replica, then concatenate in replica order
    proj_per_rep = [model.transform(f) for f in features_per_rep[site]]
    proj_concat = np.concatenate(proj_per_rep)
    projections[site] = proj_concat

    eigenvalues[site] = model.singular_values[:dim]
    timescales_ps[site] = model.timescales(lagtime=args.lag)[:dim]

    # Loadings: feature-TIC Pearson correlations on the full joint data
    feat_concat = np.concatenate(features_per_rep[site])
    loadings_mat = np.zeros((n_feat, dim))
    for j in range(n_feat):
        for i in range(dim):
            loadings_mat[j, i] = np.corrcoef(feat_concat[:, j], proj_concat[:, i])[0, 1]
    loadings[site] = loadings_mat

    print(f"== {site} ==")
    print(f"  Eigenvalues:  {np.array2string(eigenvalues[site], precision=4)}")
    print(f"  Timescales:   {np.array2string(timescales_ps[site], precision=1)} ps")
    print()


# Save TIC projections with replica labels
all_times = np.concatenate(times_per_rep)
all_replicas = np.concatenate([np.full(len(t), i + 1) for i, t in enumerate(times_per_rep)])

proj_df = pd.DataFrame({"replica": all_replicas, "time_ns": all_times})
for site in SITES:
    for i in range(dim):
        proj_df[f"{site}|TIC{i+1}"] = projections[site][:, i]
proj_df.to_csv(f"{args.prefix}_tics.csv", index=False)
print(f"Wrote {args.prefix}_tics.csv")


# Save eigenvalues and implied timescales (joint fit -> single set per site)
eig_df = pd.DataFrame({"TIC": [f"TIC{i+1}" for i in range(dim)]})
for site in SITES:
    eig_df[f"{site}|eigenvalue"] = eigenvalues[site]
    eig_df[f"{site}|timescale_ps"] = timescales_ps[site]
eig_df.to_csv(f"{args.prefix}_eigenvalues.csv", index=False)
print(f"Wrote {args.prefix}_eigenvalues.csv")


# Save TIC loadings (feature-TIC correlations)
feature_labels = [label for _, _, label in SITES["site1"]]
for dl in site_dihed_labels["site1"]:
    feature_labels.append(f"{dl}_sin")
    feature_labels.append(f"{dl}_cos")
load_df = pd.DataFrame({"feature": feature_labels})
for site in SITES:
    for i in range(dim):
        load_df[f"{site}|TIC{i+1}"] = loadings[site][:, i]
load_df.to_csv(f"{args.prefix}_loadings.csv", index=False)
print(f"Wrote {args.prefix}_loadings.csv")


# Save auxiliary distances (tracked alongside but not used in TICA)
aux_df = pd.DataFrame({"replica": all_replicas, "time_ns": all_times})
for site in aux_site_order:
    label = aux_indices[site][2]
    aux_df[f"{site}|{label}"] = np.concatenate(aux_per_rep[site])
aux_df.to_csv(f"{args.prefix}_aux.csv", index=False)
print(f"Wrote {args.prefix}_aux.csv")