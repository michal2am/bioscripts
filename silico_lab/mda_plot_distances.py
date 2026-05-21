#!/usr/bin/env python3
"""Plot atom-pair distances from a CSV produced by distances_calc.py.

Layout: rows = groups, columns = replicas. For each group, the time-series
and histogram of every replica sit in the same row so they can be compared
directly at a glance."""

import argparse
import csv
from collections import defaultdict
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="Plot atom-pair distances from CSV")
p.add_argument("-i", "--input", required=True, help="Input CSV (from distances_calc.py)")
p.add_argument("-w", "--window", type=int, default=50, help="Smoothing window in frames (default: 50, 0=off)")
p.add_argument("-b", "--bins", type=int, default=60, help="Histogram bin count (default: 60)")
p.add_argument("-o", "--prefix", default="distances", help="Output file prefix")
args = p.parse_args()

# ── Read CSV ─────────────────────────────────────────────────────────
# Header: time_ns, replica, "group|label", "group|label", ...
with open(args.input) as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = list(reader)

columns = header[2:]
groups = [c.split("|", 1)[0] for c in columns]
labels = [c.split("|", 1)[1] for c in columns]

# Group rows by replica, preserving insertion order
rep_rows = defaultdict(list)
for row in rows:
    rep_rows[row[1]].append(row)

unique_replicas = list(rep_rows.keys())

# Per-replica numpy arrays
replica_data = {}
for rep, rrows in rep_rows.items():
    times = np.array([float(r[0]) for r in rrows])
    dists = np.array([[float(x) for x in r[2:]] for r in rrows]).T  # (n_pairs, n_frames)
    replica_data[rep] = (times, dists)

# Unique groups in order of first appearance → one subplot row each
unique_groups = []
for g in groups:
    if g not in unique_groups:
        unique_groups.append(g)
group_to_row = {g: i + 1 for i, g in enumerate(unique_groups)}

n_rows = len(unique_groups)
n_cols = len(unique_replicas)

print(f"Loaded {len(labels)} pairs, {n_rows} groups, {n_cols} replicas: {unique_replicas}")
for rep, (t, _) in replica_data.items():
    print(f"  {rep}: {len(t)} frames, {t[0]:.2f}–{t[-1]:.2f} ns")


# ── Time-series figure ──────────────────────────────────────────────
def smooth(y, window=50):
    """Running average with edge handling."""
    kernel = np.ones(window) / window
    return np.convolve(y, kernel, mode="same")


fig = make_subplots(
    rows=n_rows, cols=n_cols,
    shared_xaxes="columns",  # all groups within a replica share its time axis
    #shared_yaxes="rows",     # same group compared across replicas → same y scale
    shared_yaxes="all",
    vertical_spacing=0.06, horizontal_spacing=0.04,
    row_titles=unique_groups, column_titles=unique_replicas,
)

colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

for j, (label, group) in enumerate(zip(labels, groups)):
    row_idx = group_to_row[group]
    color = colors[j % len(colors)]
    for col_idx, rep in enumerate(unique_replicas, start=1):
        times, all_dists = replica_data[rep]
        ref_d = all_dists[j, 0]
        # TODO: toggle initial value substraction
        raw = all_dists[j] - ref_d
        # raw = all_dists[j]
        show_in_legend = (col_idx == 1)  # one legend entry per pair, in first replica column

        # Raw trace (thin, semi-transparent)
        fig.add_trace(
            go.Scatter(
                x=times, y=raw, mode="lines", name=label,
                line=dict(width=0.5, color=color), opacity=0.3,
                showlegend=show_in_legend,
                legendgroup=label,
                hovertemplate=f"{label} ({rep})<br>t=%{{x:.2f}} ns<br>d=%{{y:.2f}} Å<extra></extra>",
            ),
            row=row_idx, col=col_idx,
        )
        # Smoothed trace (bold, on top)
        if args.window > 1:
            fig.add_trace(
                go.Scatter(
                    x=times, y=smooth(raw, args.window), mode="lines",
                    name=f"{label} (avg {args.window}f)", showlegend=False,
                    line=dict(width=2, color=color),
                    legendgroup=label,
                    hovertemplate=f"{label} ({rep}) smoothed<br>t=%{{x:.2f}} ns<br>d=%{{y:.2f}} Å<extra></extra>",
                ),
                row=row_idx, col=col_idx,
            )

# y-axis title only on leftmost column; x-axis title only on bottom row
for r in range(1, n_rows + 1):
    fig.update_yaxes(title_text="Distance (Å)", row=r, col=1)
for c in range(1, n_cols + 1):
    fig.update_xaxes(title_text="Time (ns)", row=n_rows, col=c)

fig.update_layout(
    height=max(900, 250 * n_rows),
    width=max(900, 450 * n_cols) + 280,  # extra width to host the right-side legend
    template="plotly_white",
    legend=dict(orientation="v", yanchor="top", y=1.0, xanchor="left", x=1.02,
                font=dict(size=9)),
    title="Atom-Pair Distances",
)

out_html = f"{args.prefix}_distances_plots.html"
fig.write_html(out_html, include_plotlyjs="cdn")
print(f"Saved {out_html}")
fig.show()

# ── Histogram figure ────────────────────────────────────────────────
# Distributions of the absolute distance for each pair, overlaid within each
# (group, replica) cell. Shared x within rows so the distance range is
# directly comparable across replicas of the same group.
fig_hist = make_subplots(
    rows=n_rows, cols=n_cols,
    shared_xaxes="rows", shared_yaxes="rows",
    vertical_spacing=0.08, horizontal_spacing=0.04,
    row_titles=unique_groups, column_titles=unique_replicas,
)

for j, (label, group) in enumerate(zip(labels, groups)):
    row_idx = group_to_row[group]
    color = colors[j % len(colors)]
    for col_idx, rep in enumerate(unique_replicas, start=1):
        _, all_dists = replica_data[rep]
        show_in_legend = (col_idx == 1)
        fig_hist.add_trace(
            go.Histogram(
                x=all_dists[j], name=label,
                marker=dict(color=color), opacity=0.5,
                nbinsx=args.bins,
                histnorm="probability density",
                showlegend=show_in_legend,
                legendgroup=label,
                hovertemplate=f"{label} ({rep})<br>d=%{{x:.2f}} Å<br>density=%{{y:.3f}}<extra></extra>",
            ),
            row=row_idx, col=col_idx,
        )

for r in range(1, n_rows + 1):
    fig_hist.update_yaxes(title_text="Density", row=r, col=1)
for c in range(1, n_cols + 1):
    fig_hist.update_xaxes(title_text="Distance (Å)", row=n_rows, col=c)

fig_hist.update_layout(
    barmode="overlay",
    height=max(900, 250 * n_rows),
    width=max(900, 450 * n_cols) + 280,
    template="plotly_white",
    legend=dict(orientation="v", yanchor="top", y=1.0, xanchor="left", x=1.02,
                font=dict(size=9)),
    title="Atom-Pair Distance Distributions",
)

out_html_hist = f"{args.prefix}_distances_histograms.html"
fig_hist.write_html(out_html_hist, include_plotlyjs="cdn")
print(f"Saved {out_html_hist}")
fig_hist.show()