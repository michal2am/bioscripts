


#!/usr/bin/env python3
"""Radial pair distribution g(r) for selected reference atoms against all
protein heavy atoms (no carbons, no hydrogens), across one or more replicas.

For each reference atom in ATOMS, distances to every target atom are histogrammed
per frame, time-averaged, normalised by spherical-shell volume, and divided by
the bulk target density (N_targets / <V_box>). Atoms belonging to the reference's
own residue are excluded from the target set so the profile reflects inter-residue
contacts only.

Layout of the plot: rows = binding sites, columns = replicas. Multiple reference
atoms sharing a binding-site label are overlaid in the same panel.
"""

import argparse
import csv
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ── Reference atom definitions ──────────────────────────────────────
# Each entry: (selection, residue_label, binding_site_label)
ATOMS = [
    ("chainID A and resid 195 and name OG1", "A:Thr195-OG1", "BS_1"),
    ("chainID C and resid 195 and name OG1", "C:Thr195-OG1", "BS_2"),

    ("chainID A and resid 196 and name NZ", "A:Lys196-NZ", "BS_1"),
    ("chainID C and resid 196 and name NZ", "C:Lys196-NZ", "BS_2"),

    ("chainID A and resid 197 and name NZ", "A:Lys197-NZ", "BS_1"),
    ("chainID C and resid 197 and name NZ", "C:Lys197-NZ", "BS_2"),

    ("chainID A and resid 201 and name OG", "A:Ser201-OG", "BS_1"),
    ("chainID C and resid 201 and name OG", "C:Ser201-OG", "BS_2"),

    ("chainID A and resid 202 and name OG1", "A:Thr202-OG1", "BS_1"),
    ("chainID C and resid 202 and name OG1", "C:THR202-OG1", "BS_2"),

    ("chainID A and resid 204 and name OG", "A:Ser204-OG", "BS_1"),
    ("chainID C and resid 204 and name OG", "C:Ser204-OG", "BS_2"),

    ("chainID A and resid 205 and name OH", "A:Tyr205-OH", "BS_1"),
    ("chainID C and resid 205 and name OH", "C:Tyr205-OH", "BS_2"),
]

# Glycan residue names to include in the target set. Common N-glycan codes
# across GLYCAM / CHARMM-CARB / GROMOS — edit to match your force field's
# naming convention if you see warnings about missing residues.
GLYCAN_RESNAMES = ["NAG", "BMA", "MAN", "AMA", "BGC", "FUC", "AFUC",
                   "SIA", "NANA", "NGA", "GAL", "BGAL", "AGAL",
                   "GLC", "AGLC", "BGLC"]

# Target: protein + glycan heavy atoms except carbons and hydrogens
TARGET_SELECTION = (
    f"(protein or resname {' '.join(GLYCAN_RESNAMES)}) "
    f"and not name C* and not name H*"
)

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="Radial pair distribution from GROMACS trajectories (multi-replica)")
p.add_argument("-s", "--tpr", required=True, help="Topology file (shared across replicas)")
p.add_argument("-f", "--xtc", required=True, nargs="+",
               help="Trajectory file(s); one per replica (labelled replica1, replica2, …)")
p.add_argument("-b", "--begin", type=float, default=None, nargs="+",
               help="Start time(s) in ns: one value for all replicas, or one per replica")
p.add_argument("-e", "--end", type=float, default=None, nargs="+",
               help="End time(s) in ns: one for all, or one per replica")
p.add_argument("-dt", "--step", type=int, default=None, nargs="+",
               help="Frame step(s): one for all, or one per replica")
p.add_argument("--rmax", type=float, default=12.0, help="Maximum distance in Å (default: 12)")
p.add_argument("--bins", type=int, default=120, help="Number of histogram bins (default: 120)")
p.add_argument("--cutoff", type=float, default=5.0,
               help="Contact cutoff in Å for the per-atom frequency plot (default: 5)")
p.add_argument("-w", "--window", type=float, default=100.0,
               help="Time window size in ns for stacked contact-frequency bars (default: 100)")
p.add_argument("-o", "--prefix", default="rdf", help="Output file prefix")
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

# Bin edges, centres, shell volumes
bin_edges = np.linspace(0.0, args.rmax, args.bins + 1)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
shell_volumes = (4.0 / 3.0) * np.pi * (bin_edges[1:] ** 3 - bin_edges[:-1] ** 3)

# ── Validate selections (once, against first replica's topology) ────
u0 = mda.Universe(args.tpr, args.xtc[0])
print(f"Topology: {u0.atoms.n_atoms} atoms; {n_replicas} replica(s)")

target_ag = u0.select_atoms(TARGET_SELECTION)
print(f"Target atoms ({TARGET_SELECTION}): {target_ag.n_atoms}")
target_idx = np.array(target_ag.indices)
target_residue_idx = np.array([u0.atoms[i].residue.resindex for i in target_idx])

# Resolve each reference atom + per-ref mask excluding same-residue targets
valid_atoms = []
ref_idx = []
ref_masks = []
for sel, res_label, bs_label in ATOMS:
    ag = u0.select_atoms(sel)
    if ag.n_atoms == 0:
        print(f"  SKIP  {res_label}  (empty selection)")
        continue
    if ag.n_atoms > 1:
        print(f"  WARN  {res_label}  using first atom (selection found {ag.n_atoms})")
    ref_atom = ag[0]
    valid_atoms.append((sel, res_label, bs_label))
    ref_idx.append(ref_atom.index)
    mask = target_residue_idx != ref_atom.residue.resindex
    ref_masks.append(mask)
    print(f"  OK    [{bs_label}]  {res_label}  → atom index {ref_atom.index}, "
          f"{mask.sum()} targets after same-residue exclusion")

ref_idx = np.array(ref_idx)

# ── Process each replica ────────────────────────────────────────────
# Per (rep, ref): accumulated histogram, frame count, summed box volume
results = []  # (rep_label, res_label, bs_label, hist, n_frames, avg_volume, n_targets)
contact_data = {}  # (rep_label, res_label) -> (contact_counts_per_target, n_frames)

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

    n_frames = len(traj_slice)
    print(f"  Total {u.trajectory.n_frames} frames; analysing {n_frames}"
          f" ({traj_slice[0].time / 1000:.2f}–{traj_slice[-1].time / 1000:.2f} ns, step={step or 1})")

    hists = [np.zeros(args.bins) for _ in valid_atoms]

    # Time windows: 0-indexed absolute, so "0-100 ns" is always the same window
    # regardless of which replica it appears in
    t_min_ns = traj_slice[0].time / 1000
    t_max_ns = traj_slice[-1].time / 1000
    first_win = int(np.floor(t_min_ns / args.window))
    last_win  = int(np.floor(t_max_ns / args.window))
    n_wins = last_win - first_win + 1
    window_labels = [f"{(first_win + w) * args.window:.0f}-{(first_win + w + 1) * args.window:.0f} ns"
                     for w in range(n_wins)]
    window_contact_counts = [np.zeros((n_wins, int(m.sum())), dtype=int) for m in ref_masks]
    sum_volume = 0.0

    for i, ts in enumerate(traj_slice):
        pos = ts.positions
        win_idx = max(0, min(n_wins - 1, int(ts.time / 1000 / args.window) - first_win))
        # PBC-aware distances: (n_refs × n_targets)
        d = distance_array(pos[ref_idx], pos[target_idx], box=ts.dimensions)
        for k in range(len(ref_idx)):
            dists = d[k][ref_masks[k]]  # exclude same residue
            counts, _ = np.histogram(dists, bins=bin_edges)
            hists[k] += counts
            window_contact_counts[k][win_idx] += (dists < args.cutoff).astype(int)
        sum_volume += ts.volume
        if (i + 1) % 2000 == 0 or i == n_frames - 1:
            print(f"  Frame {i + 1}/{n_frames}")

    avg_volume = sum_volume / n_frames if n_frames > 0 else np.nan
    for k, (sel, res_label, bs_label) in enumerate(valid_atoms):
        results.append((rep_label, res_label, bs_label,
                        hists[k], n_frames, avg_volume, int(ref_masks[k].sum())))
        contact_data[(rep_label, res_label)] = (window_contact_counts[k], n_frames, window_labels)

# ── Compute g(r) ────────────────────────────────────────────────────
# g(r) = (counts/n_frames) / (shell_volume * N_targets / <V_box>)
print(f"\nComputing g(r)...")
rdf_data = {}  # (rep, res_label) -> g_r
for rep, res_label, bs_label, hist, n_frames, avg_volume, n_targets in results:
    rho_bulk = n_targets / avg_volume
    g_r = (hist / n_frames) / (shell_volumes * rho_bulk)
    rdf_data[(rep, res_label)] = g_r

# ── Save CSV (long format) ──────────────────────────────────────────
# Columns: r_ang, replica, "bs|res" for each reference
unique_replicas = [f"replica{i+1}" for i in range(n_replicas)]
header = ["r_ang", "replica"] + [f"{bs}|{res}" for _, res, bs in valid_atoms]
out_csv = f"{args.prefix}_rdf.csv"
with open(out_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(header)
    for rep in unique_replicas:
        for i, r in enumerate(bin_centers):
            row = [f"{r:.4f}", rep]
            for _, res, _ in valid_atoms:
                if (rep, res) in rdf_data:
                    row.append(f"{rdf_data[(rep, res)][i]:.6f}")
                else:
                    row.append("nan")
            w.writerow(row)
print(f"Saved {out_csv}")

# ── Plot ─────────────────────────────────────────────────────────────
# Layout: rows = binding sites, cols = replicas
unique_bs = []
for _, _, bs in valid_atoms:
    if bs not in unique_bs:
        unique_bs.append(bs)
bs_to_row = {bs: i + 1 for i, bs in enumerate(unique_bs)}

n_rows = len(unique_bs)
n_cols = len(unique_replicas) + 1  # +1 for average column

fig = make_subplots(
    rows=n_rows, cols=n_cols,
    shared_xaxes="all", shared_yaxes="rows",
    vertical_spacing=0.04, horizontal_spacing=0.01,
    row_titles=unique_bs, column_titles=list(unique_replicas) + ["Average"],
)

colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

for j, (sel, res_label, bs_label) in enumerate(valid_atoms):
    row_idx = bs_to_row[bs_label]
    color = colors[j % len(colors)]
    for c_idx, rep in enumerate(unique_replicas, start=1):
        if (rep, res_label) not in rdf_data:
            continue
        g_r = rdf_data[(rep, res_label)]
        show_in_legend = (c_idx == 1)
        fig.add_trace(
            go.Scatter(
                x=bin_centers, y=g_r, mode="lines",
                line=dict(width=2, color=color),
                name=res_label,
                showlegend=show_in_legend,
                legendgroup=res_label,
                hovertemplate=f"{res_label} ({rep})<br>r=%{{x:.2f}} Å<br>g(r)=%{{y:.2f}}<extra></extra>",
            ),
            row=row_idx, col=c_idx,
        )

    # Average g(r) across replicas in the last column
    g_rs = [rdf_data[(rep, res_label)] for rep in unique_replicas
            if (rep, res_label) in rdf_data]
    if g_rs:
        g_r_avg = np.mean(g_rs, axis=0)
        fig.add_trace(
            go.Scatter(
                x=bin_centers, y=g_r_avg, mode="lines",
                line=dict(width=2, color=color),
                name=res_label,
                showlegend=False,
                legendgroup=res_label,
                hovertemplate=f"{res_label} (avg)<br>r=%{{x:.2f}} Å<br>g(r)=%{{y:.2f}}<extra></extra>",
            ),
            row=row_idx, col=n_cols,
        )

# Reference line at g(r) = 1 (bulk density) on every panel
for r in range(1, n_rows + 1):
    for c in range(1, n_cols + 1):
        fig.add_hline(y=1.0, line=dict(width=1, color="gray", dash="dot"),
                      row=r, col=c)

# Y-axis title on leftmost column; x-axis title on bottom row
for r in range(1, n_rows + 1):
    fig.update_yaxes(title_text="g(r)", row=r, col=1)
for c in range(1, n_cols + 1):
    fig.update_xaxes(title_text="r (Å)", row=n_rows, col=c)

fig.update_layout(
    height=max(600, 280 * n_rows),
    width=max(900, 300 * n_cols) + 280,
    template="plotly_white",
    legend=dict(orientation="v", yanchor="top", y=1.0, xanchor="left", x=1.02,
                font=dict(size=9)),
    title=f"Radial Pair Distribution g(r) — vs {TARGET_SELECTION}",
)

out_html = f"{args.prefix}_rdf_plots.html"
fig.write_html(out_html, include_plotlyjs="cdn")
print(f"Saved {out_html}")
fig.show()

# ── Contact-frequency figure ────────────────────────────────────────
# For each reference atom: bar chart of contact frequency per target atom.
# Targets are filtered to those that came within --cutoff Å in ≥1 replica,
# and shown in every replica column for that ref (so 0-bars highlight
# replicas where a contact partner was absent).
print(f"\nBuilding contact-frequency plot (cutoff {args.cutoff} Å)...")

# Per-ref absolute target indices (after same-residue exclusion)
ref_target_indices = [target_idx[m] for m in ref_masks]

per_ref_selected = {}  # res_label -> (sel_local indices, target_labels)
for k, (sel, res_label, bs_label) in enumerate(valid_atoms):
    n_t = len(ref_target_indices[k])
    total_freq = np.zeros(n_t)
    any_contact = np.zeros(n_t, dtype=bool)
    for rep in unique_replicas:
        if (rep, res_label) in contact_data:
            window_counts, n_frames, _ = contact_data[(rep, res_label)]
            counts = window_counts.sum(axis=0)
            freq = counts / n_frames
            total_freq += freq
            any_contact |= (freq > 0)
    sel_local = np.where(any_contact)[0]
    # Sort by total frequency across replicas, descending
    sel_local = sel_local[np.argsort(-total_freq[sel_local])]
    target_labels = []
    for idx_local in sel_local:
        a = u0.atoms[ref_target_indices[k][idx_local]]
        target_labels.append(f"{a.chainID}:{a.resname}{a.resid}-{a.name}")
    per_ref_selected[res_label] = (sel_local, target_labels)
    print(f"  {res_label}: {len(sel_local)} target atom(s) within {args.cutoff} Å in ≥1 replica")

n_rows_c = len(valid_atoms)
n_cols_c = len(unique_replicas) + 1  # +1 for average column

# Map each window label to a consistent colour across replicas (by absolute time)
from plotly.colors import sample_colorscale
all_window_labels = set()
for (_, _), (_, _, win_lbls) in contact_data.items():
    all_window_labels.update(win_lbls)
sorted_window_labels = sorted(all_window_labels, key=lambda s: float(s.split("-")[0]))
n_total_wins = len(sorted_window_labels)
window_colors = {
    lbl: sample_colorscale("Viridis", [i / max(n_total_wins - 1, 1)])[0]
    for i, lbl in enumerate(sorted_window_labels)
}

fig_contacts = make_subplots(
    rows=n_rows_c, cols=n_cols_c,
    shared_xaxes="rows",   # same target set across replicas of a given ref
    shared_yaxes="all",    # frequency ∈ [0, 1]
    vertical_spacing=0.01, horizontal_spacing=0.01,
    row_titles=[res for _, res, _ in valid_atoms],
    column_titles=list(unique_replicas) + ["Average"],
)

shown_windows = set()  # track which window labels have appeared in legend
for k, (sel, res_label, bs_label) in enumerate(valid_atoms):
    row = k + 1
    sel_local, target_labels = per_ref_selected[res_label]
    if len(sel_local) == 0:
        continue
    for c_idx, rep in enumerate(unique_replicas, start=1):
        if (rep, res_label) not in contact_data:
            continue
        window_counts, n_frames, win_lbls = contact_data[(rep, res_label)]
        for w_idx, win_lbl in enumerate(win_lbls):
            seg = window_counts[w_idx][sel_local] / n_frames
            show_in_legend = win_lbl not in shown_windows
            shown_windows.add(win_lbl)
            fig_contacts.add_trace(
                go.Bar(
                    x=target_labels, y=seg,
                    marker=dict(color=window_colors[win_lbl]),
                    name=win_lbl,
                    showlegend=show_in_legend,
                    legendgroup=win_lbl,
                    hovertemplate=(f"{res_label} ({rep})<br>"
                                   f"partner: %{{x}}<br>"
                                   f"window: {win_lbl}<br>"
                                   f"segment = %{{y:.3f}}<extra></extra>"),
                ),
                row=row, col=c_idx,
            )

    # Average column: per window, mean segment across replicas that include it
    for win_lbl in sorted_window_labels:
        seg_heights = []
        for rep in unique_replicas:
            if (rep, res_label) not in contact_data:
                continue
            window_counts, n_frames, win_lbls = contact_data[(rep, res_label)]
            if win_lbl in win_lbls:
                w_idx = win_lbls.index(win_lbl)
                seg_heights.append(window_counts[w_idx][sel_local] / n_frames)
        if not seg_heights:
            continue
        avg_seg = np.mean(seg_heights, axis=0)
        fig_contacts.add_trace(
            go.Bar(
                x=target_labels, y=avg_seg,
                marker=dict(color=window_colors[win_lbl]),
                name=win_lbl,
                showlegend=False,
                legendgroup=win_lbl,
                hovertemplate=(f"{res_label} (avg)<br>"
                               f"partner: %{{x}}<br>"
                               f"window: {win_lbl}<br>"
                               f"segment = %{{y:.3f}}<extra></extra>"),
            ),
            row=row, col=n_cols_c,
        )

# Y-axis: frequency in [0, 1.05], title on leftmost column
for r in range(1, n_rows_c + 1):
    fig_contacts.update_yaxes(title_text="Frequency", range=[0, 1.05], row=r, col=1)

# X-axis: tilt long target labels on every panel
for r in range(1, n_rows_c + 1):
    for c in range(1, n_cols_c + 1):
        fig_contacts.update_xaxes(tickangle=-45, row=r, col=c)

fig_contacts.update_layout(
    barmode="stack",
    height=max(600, 350 * n_rows_c),
    width=max(900, 300 * n_cols_c) + 280,
    template="plotly_white",
    legend=dict(orientation="v", yanchor="top", y=1.0, xanchor="left", x=1.02,
                font=dict(size=9), title=dict(text="Time window")),
    title=f"Contact frequencies within {args.cutoff} Å (stacked by {args.window:.0f} ns window)",
)

out_html_contacts = f"{args.prefix}_contacts.html"
fig_contacts.write_html(out_html_contacts, include_plotlyjs="cdn")
print(f"Saved {out_html_contacts}")
fig_contacts.show()