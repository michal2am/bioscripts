#!/usr/bin/env python3
"""
TICA on B9-B10 antiparallel beta strand H-bond distances.
Separate TICA per binding site (chain A = site1, chain C = site2).

Requires: numpy, pandas, MDAnalysis, deeptime
"""

import argparse
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_bonds
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
        ("chainID A and resid 200 and name CG", "chainID B and resid 46 and name CG", "200-CG↔46-CG"),
        # ("chainID C and resid 200 and name CG", "chainID D and resid 46 and name CG", "NEI200-CG↔46-CG")

    ],
    "site2": [
        # GABA
        ("chainID C and resid 202 and name OG1", "chainID I and resname ABU and name C4", "Thr202-OG1↔ABU-C4"),
        ("chainID D and resid 67 and name CZ", "chainID I and resname ABU and name C4", "Arg67-CZ↔ABU-C4"),
        ("chainID D and resid 130 and name OG1", "chainID I and resname ABU and name C4", "Thr130-OG1↔ABU-C4"),
        ("chainID C and resid 155 and name CD", "chainID I and resname ABU and name N", "Glu155-CD↔ABU-N"),
        ("chainID C and resid 97 and name OH", "chainID I and resname ABU and name N", "Tyr97-OH↔ABU-N"),
        ("chainID C and resid 156 and name O", "chainID I and resname ABU and name N", "Ser156-O↔ABU-N"),
        # h-bonds
        ("chainID C and resid 207 and name N", "chainID C and resid 196 and name O", "207-N↔196-O"),
        ("chainID C and resid 207 and name O", "chainID C and resid 196 and name N", "207-O↔196-N"),
        ("chainID C and resid 205 and name N", "chainID C and resid 198 and name O", "205-N↔198-O"),
        ("chainID C and resid 205 and name O", "chainID C and resid 198 and name N", "205-O↔198-N"),
        ("chainID C and resid 203 and name N", "chainID C and resid 200 and name O", "203-N↔200-O"),
        ("chainID C and resid 203 and name O", "chainID C and resid 200 and name N", "203-O↔200-N"),
        # trio
        ("chainID C and resid 153 and name CD", "chainID C and resid 207 and name NH1", "153-CD↔207-NH1"),
        ("chainID C and resid 155 and name CD", "chainID C and resid 207 and name NH2", "155-CD↔207-NH2"),
        # cap h-bond
        ("chainID C and resid 205 and name OH", "chainID C and resid 202 and name OG1", "205-OH↔202-OG1"),
        # interfaces
        ("chainID C and resid 200 and name CG", "chainID D and resid 46 and name CG", "200-CG↔46-CG"),
        # ("chainID A and resid 200 and name CG", "chainID B and resid 46 and name CG", "NEI200-CG↔46-CG")

    ],
}


p = argparse.ArgumentParser(description="TICA on B9-B10 H-bond distances, per binding site")
p.add_argument("-s", "--tpr", required=True, help="Topology file (.tpr/.gro/.pdb)")
p.add_argument("-f", "--xtc", required=True, help="Trajectory file (.xtc/.trr)")
p.add_argument("-b", "--begin", type=float, default=None, help="Start time in ns (default: first frame)")
p.add_argument("-e", "--end", type=float, default=None, help="End time in ns (default: last frame)")
p.add_argument("-dt", "--step", type=int, default=None, help="Process every Nth frame (default: all)")
p.add_argument("-lag", "--lag", type=float, default=1000.0, help="TICA lag time in ps (default: 1000)")
p.add_argument("-n", "--n-tics", type=int, default=None, help="Number of TICs to keep (default: all features)")
p.add_argument("-o", "--prefix", default="tica", help="Output file prefix")
args = p.parse_args()


u = mda.Universe(args.tpr, args.xtc)
dt_ps = u.trajectory.dt
step = args.step or 1
effective_dt_ps = dt_ps * step
lag_frames = max(1, int(round(args.lag / effective_dt_ps)))

# Resolve atom indices once
site_indices = {}
for site, pairs in SITES.items():
    idx1 = [u.select_atoms(s1).indices[0] for s1, _, _ in pairs]
    idx2 = [u.select_atoms(s2).indices[0] for _, s2, _ in pairs]
    site_indices[site] = (np.array(idx1), np.array(idx2))

# Collect distances
times_ns = []
features = {site: [] for site in SITES}

for ts in u.trajectory:
    t_ns = ts.time / 1000.0
    if args.begin is not None and t_ns < args.begin:
        continue
    if args.end is not None and t_ns > args.end:
        break
    if (ts.frame % step) != 0:
        continue
    times_ns.append(t_ns)
    positions = u.atoms.positions
    for site, (i1, i2) in site_indices.items():
        d = calc_bonds(positions[i1], positions[i2], box=ts.dimensions)
        features[site].append(d)

times_ns = np.array(times_ns)
for site in features:
    features[site] = np.array(features[site])

print(f"Loaded {len(times_ns)} frames (effective dt = {effective_dt_ps} ps)")
print(f"Lag time: {args.lag} ps = {lag_frames} frames")
print()

# Fit TICA per site
n_feat = len(SITES["site1"])
dim = args.n_tics if args.n_tics is not None else n_feat
dim = min(dim, n_feat)

projections = {}
eigenvalues = {}
timescales_ps = {}
loadings = {}

for site in SITES:
    tica = TICA(lagtime=lag_frames, dim=dim)
    model = tica.fit(features[site]).fetch_model()
    proj = model.transform(features[site])
    projections[site] = proj
    eigenvalues[site] = model.singular_values[:dim]
    timescales_ps[site] = model.timescales(lagtime=args.lag)[:dim]
    # Loadings: feature-TIC Pearson correlations
    loadings_mat = np.zeros((n_feat, dim))
    for j in range(n_feat):
        for i in range(dim):
            loadings_mat[j, i] = np.corrcoef(features[site][:, j], proj[:, i])[0, 1]
    loadings[site] = loadings_mat

    print(f"== {site} ==")
    print(f"  Eigenvalues:  {np.array2string(eigenvalues[site], precision=4)}")
    print(f"  Timescales:   {np.array2string(timescales_ps[site], precision=1)} ps")
    print()


# Save TIC projections (time series)
proj_df = pd.DataFrame({"time_ns": times_ns})
for site in SITES:
    for i in range(dim):
        proj_df[f"{site}|TIC{i+1}"] = projections[site][:, i]
proj_df.to_csv(f"{args.prefix}_tics.csv", index=False)
print(f"Wrote {args.prefix}_tics.csv")


# Save eigenvalues and implied timescales
eig_df = pd.DataFrame({"TIC": [f"TIC{i+1}" for i in range(dim)]})
for site in SITES:
    eig_df[f"{site}|eigenvalue"] = eigenvalues[site]
    eig_df[f"{site}|timescale_ps"] = timescales_ps[site]
eig_df.to_csv(f"{args.prefix}_eigenvalues.csv", index=False)
print(f"Wrote {args.prefix}_eigenvalues.csv")


# Save TIC loadings (feature-TIC correlations)
feature_labels = [label for _, _, label in SITES["site1"]]
load_df = pd.DataFrame({"feature": feature_labels})
for site in SITES:
    for i in range(dim):
        load_df[f"{site}|TIC{i+1}"] = loadings[site][:, i]
load_df.to_csv(f"{args.prefix}_loadings.csv", index=False)
print(f"Wrote {args.prefix}_loadings.csv")