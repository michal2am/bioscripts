#!/usr/bin/env python3
"""Calculate RMSD (total + per chain) and RMSF per residue from a GROMACS trajectory."""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="RMSD & RMSF from GROMACS trajectory")
p.add_argument("-s", "--tpr", required=True, help="Topology file (.tpr/.gro/.pdb)")
p.add_argument("-f", "--xtc", required=True, help="Trajectory file (.xtc/.trr)")
p.add_argument("-sel", "--selection", default="backbone", help="Atom selection for RMSD/RMSF (default: backbone)")
p.add_argument("-ref", "--ref-frame", type=int, default=0, help="Reference frame index (default: 0)")
p.add_argument("-o", "--prefix", default="analysis", help="Output file prefix")
args = p.parse_args()

# ── Load universe ────────────────────────────────────────────────────
u = mda.Universe(args.tpr, args.xtc)
sel = args.selection
ref_frame = args.ref_frame

print(f"Loaded: {u.trajectory.n_frames} frames, {u.atoms.n_atoms} atoms")
print(f"Selection: '{sel}'  |  Reference frame: {ref_frame}")

# ── Identify chains ─────────────────────────────────────────────────
seg_ids = sorted(set(u.select_atoms(sel).segids))
if len(seg_ids) <= 1:
    # fallback: try chainIDs
    seg_ids = sorted(set(u.select_atoms(sel).chainIDs))
    chain_key = "chainID"
else:
    chain_key = "segid"

print(f"Chains found ({chain_key}): {seg_ids}")

# ── 1. Total RMSD ───────────────────────────────────────────────────
R_total = rms.RMSD(u, u, select=sel, ref_frame=ref_frame)
R_total.run(verbose=True)

time_ns = R_total.results.rmsd[:, 1] / 1000.0  # ps → ns
rmsd_total = R_total.results.rmsd[:, 2]  # Å

# ── 2. Per-chain RMSD ───────────────────────────────────────────────
chain_rmsd = {}
for cid in seg_ids:
    chain_sel = f"{sel} and {chain_key} {cid}"
    if u.select_atoms(chain_sel).n_atoms == 0:
        continue
    R = rms.RMSD(u, u, select=chain_sel, ref_frame=ref_frame)
    R.run(verbose=False)
    chain_rmsd[cid] = R.results.rmsd[:, 2]
    print(f"  Chain {cid}: mean RMSD = {chain_rmsd[cid].mean():.3f} Å")

# ── 3. RMSF per residue ─────────────────────────────────────────────
# Align trajectory to reference first
align.AlignTraj(u, u, select=sel, ref_frame=ref_frame, in_memory=True).run(verbose=True)

atoms = u.select_atoms(sel)
from MDAnalysis.analysis.rms import RMSF as calcRMSF

rmsf_calc = calcRMSF(atoms).run(verbose=True)

# Map atom RMSF → per-residue (mean over backbone atoms in each residue)
residues = atoms.residues
resids = residues.resids
resnames = residues.resnames
seg_per_res = residues.segids if chain_key == "segid" else [atoms.select_atoms(f"resid {r}").chainIDs[0] for r in
                                                            resids]

rmsf_per_res = np.array([
    rmsf_calc.results.rmsf[atoms.resindices == res.resindex].mean()
    for res in residues
])

# ── Save CSV outputs ─────────────────────────────────────────────────
# RMSD
rmsd_header = "time_ns,rmsd_total_A," + ",".join(f"rmsd_{c}_A" for c in chain_rmsd)
rmsd_data = np.column_stack([time_ns, rmsd_total] + [chain_rmsd[c] for c in chain_rmsd])
np.savetxt(f"{args.prefix}_rmsd.csv", rmsd_data, delimiter=",", header=rmsd_header, comments="")
print(f"Saved {args.prefix}_rmsd.csv")

# RMSF
with open(f"{args.prefix}_rmsf.csv", "w") as f:
    f.write("chain,resid,resname,rmsf_A\n")
    for i, res in enumerate(residues):
        f.write(f"{seg_per_res[i]},{resids[i]},{resnames[i]},{rmsf_per_res[i]:.4f}\n")
print(f"Saved {args.prefix}_rmsf.csv")


'''
# ── Plots ────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(12, 8), constrained_layout=True)

# RMSD plot
ax = axes[0]
ax.plot(time_ns, rmsd_total, "k-", lw=1.2, label="Total")
for cid, vals in chain_rmsd.items():
    ax.plot(time_ns, vals, lw=0.9, alpha=0.8, label=f"Chain {cid}")
ax.set_xlabel("Time (ns)")
ax.set_ylabel("RMSD (Å)")
ax.set_title("RMSD")
ax.legend(fontsize=8)

# RMSF plot
ax = axes[1]
unique_chains = sorted(set(seg_per_res))
offset = 0
tick_pos, tick_lab = [], []
for cid in unique_chains:
    mask = np.array(seg_per_res) == cid
    x = np.arange(mask.sum()) + offset
    ax.bar(x, rmsf_per_res[mask], width=1.0, alpha=0.75, label=f"Chain {cid}")
    # tick every 10th residue
    ids = resids[mask]
    for j in range(0, len(ids), 10):
        tick_pos.append(x[j])
        tick_lab.append(str(ids[j]))
    offset += mask.sum() + 3  # small gap between chains

ax.set_xticks(tick_pos)
ax.set_xticklabels(tick_lab, rotation=90, fontsize=6)
ax.set_xlabel("Residue ID")
ax.set_ylabel("RMSF (Å)")
ax.set_title("RMSF per residue")
ax.legend(fontsize=8)

plt.savefig(f"{args.prefix}_plots.png", dpi=200)
print(f"Saved {args.prefix}_plots.png")
plt.show()
'''

# ── Interactive Plots (Plotly) ────────────────────────────────────────
fig = make_subplots(rows=2, cols=1, subplot_titles=("RMSD", "RMSF per residue"),
                    vertical_spacing=0.12)

# RMSD traces
fig.add_trace(go.Scatter(x=time_ns, y=rmsd_total, mode="lines", name="Total",
                         line=dict(color="black", width=1.5)), row=1, col=1)
for cid, vals in chain_rmsd.items():
    fig.add_trace(go.Scatter(x=time_ns, y=vals, mode="lines", name=f"Chain {cid}",
                             line=dict(width=1.2)), row=1, col=1)

fig.update_xaxes(title_text="Time (ns)", row=1, col=1)
fig.update_yaxes(title_text="RMSD (Å)", row=1, col=1)

# RMSF bar traces
unique_chains = sorted(set(seg_per_res))
seg_arr = np.array(seg_per_res)
for cid in unique_chains:
    mask = seg_arr == cid
    labels = [f"{resnames[i]}{resids[i]}" for i in np.where(mask)[0]]
    fig.add_trace(go.Bar(x=labels, y=rmsf_per_res[mask], name=f"Chain {cid}",
                         hovertemplate="%{x}: %{y:.3f} Å<extra>Chain " + str(cid) + "</extra>"),
                  row=2, col=1)

fig.update_xaxes(title_text="Residue", tickangle=-90, tickfont=dict(size=8), row=2, col=1)
fig.update_yaxes(title_text="RMSF (Å)", row=2, col=1)

fig.update_layout(height=900, template="plotly_white", barmode="group",
                  legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))

out_html = f"{args.prefix}_rmsd_rmsf_plots.html"
fig.write_html(out_html, include_plotlyjs="cdn")
print(f"Saved {out_html}")
fig.show()