#!/usr/bin/env python3
"""Plot per-residue dihedral angles from a CSV produced by dihedrals_calc.py.

Layout: rows = residues, columns = replicas. For each residue, every replica
sits in the same row so they can be compared at a glance.

Raw data is shown as scatter markers (no transformations), with a circular-mean
smoothing line on top. The smoothed line uses NaN-break at ±180° wraps to
avoid drawing fake ~360° jumps. A second figure shows angle-value histograms.
"""

import argparse
import csv
from collections import defaultdict
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="Plot dihedral angles from CSV")
p.add_argument("-i", "--input", required=True, help="Input CSV (from dihedrals_calc.py)")
p.add_argument("-w", "--window", type=int, default=50, help="Smoothing window in frames (default: 50, 0=off)")
p.add_argument("-t", "--threshold", type=float, default=180,
               help="Wrap-detection threshold in deg for the smoothed line (default: 180)")
p.add_argument("-a", "--angles", nargs="+", default=None,
               help="Angle types to plot, e.g. 'phi psi' for backbone only (default: all)")
p.add_argument("-b", "--bins", type=int, default=60, help="Histogram bin count (default: 60)")
p.add_argument("-o", "--prefix", default="dihedrals", help="Output file prefix")
args = p.parse_args()

# ── Read CSV ─────────────────────────────────────────────────────────
# Header: time_ns, replica, "group|label", "group|label", ...
with open(args.input) as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = list(reader)

columns = header[2:]
groups_all = [c.split("|", 1)[0] for c in columns]
labels_all = [c.split("|", 1)[1] for c in columns]

# Optional angle-type filter — label looks like "A:200:phi", so the type is
# the part after the last colon.
if args.angles is not None:
    keep = [i for i, l in enumerate(labels_all) if l.rsplit(":", 1)[-1] in args.angles]
else:
    keep = list(range(len(labels_all)))

groups = [groups_all[i] for i in keep]
labels = [labels_all[i] for i in keep]

# Group rows by replica, preserving insertion order
rep_rows = defaultdict(list)
for row in rows:
    rep_rows[row[1]].append(row)

unique_replicas = list(rep_rows.keys())

# Per-replica numpy arrays (subset to the kept angle series)
replica_data = {}
for rep, rrows in rep_rows.items():
    times = np.array([float(r[0]) for r in rrows])
    angles = np.array([[float(x) for x in r[2:]] for r in rrows]).T  # (all_series, n_frames)
    replica_data[rep] = (times, angles[keep])

# Unique groups in order of first appearance → one subplot row each
unique_groups = []
for g in groups:
    if g not in unique_groups:
        unique_groups.append(g)
group_to_row = {g: i + 1 for i, g in enumerate(unique_groups)}

n_rows = len(unique_groups)
n_cols = len(unique_replicas)

print(f"Loaded {len(labels)} angles ({'all' if args.angles is None else ' '.join(args.angles)}), "
      f"{n_rows} groups, {n_cols} replicas: {unique_replicas}")
for rep, (t, _) in replica_data.items():
    print(f"  {rep}: {len(t)} frames, {t[0]:.2f}–{t[-1]:.2f} ns")


# ── Helpers ─────────────────────────────────────────────────────────
def break_wraps(t, y, threshold=180):
    """Insert NaN where |Δy| > threshold so Plotly breaks the line there."""
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)
    jumps = np.where(np.abs(np.diff(y)) > threshold)[0]
    return np.insert(t, jumps + 1, np.nan), np.insert(y, jumps + 1, np.nan)


def circular_smooth(y_deg, window):
    """Running circular mean of an angle series (in degrees)."""
    rad = np.radians(y_deg)
    kernel = np.ones(window) / window
    s = np.convolve(np.sin(rad), kernel, mode="same")
    c = np.convolve(np.cos(rad), kernel, mode="same")
    return np.degrees(np.arctan2(s, c))


colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

# ── Time-series figure ──────────────────────────────────────────────
fig = make_subplots(
    rows=n_rows, cols=n_cols,
    shared_xaxes="columns",  # groups within a replica share its time axis
    shared_yaxes="rows",     # angle range is fixed anyway, but keeps things tidy
    vertical_spacing=0.004, horizontal_spacing=0.04,
    row_titles=unique_groups, column_titles=unique_replicas,
)

for j, (label, group) in enumerate(zip(labels, groups)):
    row_idx = group_to_row[group]
    color = colors[j % len(colors)]

    if args.window > 1:
        opa = 0.3
    else:
        opa = 1

    for col_idx, rep in enumerate(unique_replicas, start=1):
        times, all_angles = replica_data[rep]
        angles = all_angles[j]
        show_in_legend = (col_idx == 1)  # one legend entry per pair, in first replica column

        # Raw scatter (markers only — no line so no wrap artefacts)
        fig.add_trace(
            go.Scatter(
                x=times, y=angles, mode="markers", name=label,
                marker=dict(size=3, color=color), opacity=opa,
                showlegend=show_in_legend,
                legendgroup=label,
                hovertemplate=f"{label} ({rep})<br>t=%{{x:.2f}} ns<br>θ=%{{y:.1f}}°<extra></extra>",
            ),
            row=row_idx, col=col_idx,
        )
        # Smoothed line via circular mean (break on wraps)
        if args.window > 1:
            smoothed = circular_smooth(angles, args.window)
            t_sm, y_sm = break_wraps(times, smoothed, args.threshold)
            fig.add_trace(
                go.Scatter(
                    x=t_sm, y=y_sm, mode="lines",
                    name=f"{label} (avg {args.window}f)", showlegend=False,
                    line=dict(width=2, color=color),
                    legendgroup=label,
                    hovertemplate=f"{label} ({rep}) smoothed<br>t=%{{x:.2f}} ns<br>θ=%{{y:.1f}}°<extra></extra>",
                ),
                row=row_idx, col=col_idx,
            )

# Y-axis: fixed range + tick marks on every row; title only on leftmost column
for r in range(1, n_rows + 1):
    fig.update_yaxes(range=[-200, 200],
                     tickvals=[-180, -120, -60, 0, 60, 120, 180],
                     row=r, col=1)
    fig.update_yaxes(title_text="Angle (°)", row=r, col=1)

# X-axis title only on bottom row
for c in range(1, n_cols + 1):
    fig.update_xaxes(title_text="Time (ns)", row=n_rows, col=c)

fig.update_layout(
    height=max(900, 250 * n_rows),
    width=max(900, 450 * n_cols) + 280,  # extra width to host the right-side legend
    template="plotly_white",
    legend=dict(orientation="v", yanchor="top", y=1.0, xanchor="left", x=1.02,
                font=dict(size=9)),
    title="Dihedral Angles",
)

out_html = f"{args.prefix}_dihedrals_plots.html"
fig.write_html(out_html, include_plotlyjs="cdn")
print(f"Saved {out_html}")
fig.show()

# ── Histogram figure ────────────────────────────────────────────────
# Angle-value distributions, overlaid within each (residue, replica) cell.
# Bins are pinned to [-180, 180] so every panel and replica shares the same
# binning and the distributions are directly comparable.
fig_hist = make_subplots(
    rows=n_rows, cols=n_cols,
    shared_xaxes="all", shared_yaxes="rows",
    vertical_spacing=0.004, horizontal_spacing=0.04,
    row_titles=unique_groups, column_titles=unique_replicas,
)

for j, (label, group) in enumerate(zip(labels, groups)):
    row_idx = group_to_row[group]
    color = colors[j % len(colors)]
    for col_idx, rep in enumerate(unique_replicas, start=1):
        _, all_angles = replica_data[rep]
        show_in_legend = (col_idx == 1)
        fig_hist.add_trace(
            go.Histogram(
                x=all_angles[j], name=label,
                marker=dict(color=color), opacity=0.5,
                xbins=dict(start=-180, end=180, size=360 / args.bins),
                histnorm="probability density",
                showlegend=show_in_legend,
                legendgroup=label,
                hovertemplate=f"{label} ({rep})<br>θ=%{{x:.1f}}°<br>density=%{{y:.4f}}<extra></extra>",
            ),
            row=row_idx, col=col_idx,
        )

# Y-axis title on leftmost column; x-axis fixed to [-180, 180] with rotamer ticks
for r in range(1, n_rows + 1):
    fig_hist.update_yaxes(title_text="Density", row=r, col=1)
    for c in range(1, n_cols + 1):
        fig_hist.update_xaxes(range=[-180, 180],
                              tickvals=[-180, -120, -60, 0, 60, 120, 180],
                              row=r, col=c)
for c in range(1, n_cols + 1):
    fig_hist.update_xaxes(title_text="Angle (°)", row=n_rows, col=c)

fig_hist.update_layout(
    barmode="overlay",
    height=max(900, 250 * n_rows),
    width=max(900, 450 * n_cols) + 280,
    template="plotly_white",
    legend=dict(orientation="v", yanchor="top", y=1.0, xanchor="left", x=1.02,
                font=dict(size=9)),
    title="Dihedral Angle Distributions",
)

out_html_hist = f"{args.prefix}_dihedrals_histograms.html"
fig_hist.write_html(out_html_hist, include_plotlyjs="cdn")
print(f"Saved {out_html_hist}")
fig_hist.show()