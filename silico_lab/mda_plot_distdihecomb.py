#!/usr/bin/env python3
"""Combined grid plot of distances and dihedrals.

Reads two long-format CSVs (from distances_calc.py and dihedrals_calc.py),
stacks distance groups on top of dihedral groups as separate rows, with
replicas as columns and shared time within each column.

Distances: thin raw line + smoothed bold line (Δ from initial frame).
Dihedrals: scatter markers + circular-mean smoothed line, NaN-broken at ±180°.
"""

import argparse
import csv
from collections import defaultdict
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="Combined distance + dihedral plot from two CSVs")
p.add_argument("-d", "--distances", required=True, help="Distance CSV (from distances_calc.py)")
p.add_argument("-a", "--angles",    required=True, help="Dihedral CSV (from dihedrals_calc.py)")
p.add_argument("-w", "--window",    type=int,   default=50,  help="Smoothing window in frames (default: 50, 0=off)")
p.add_argument("-t", "--threshold", type=float, default=180,
               help="Wrap threshold for dihedral smoothed line (default: 180)")
p.add_argument("-o", "--prefix", default="combined", help="Output file prefix")
args = p.parse_args()


# ── CSV reader ──────────────────────────────────────────────────────
def read_long_csv(path):
    """Read a long-format CSV with header: time_ns, replica, "group|label", ...
    Returns (groups, labels, replica_data) where replica_data is
    {rep_label: (times, values)} and values is (n_series, n_frames)."""
    with open(path) as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)
    columns = header[2:]
    groups = [c.split("|", 1)[0] for c in columns]
    labels = [c.split("|", 1)[1] for c in columns]
    rep_rows = defaultdict(list)
    for row in rows:
        rep_rows[row[1]].append(row)
    replica_data = {}
    for rep, rrows in rep_rows.items():
        times = np.array([float(r[0]) for r in rrows])
        vals  = np.array([[float(x) for x in r[2:]] for r in rrows]).T
        replica_data[rep] = (times, vals)
    return groups, labels, replica_data


d_groups, d_labels, d_data = read_long_csv(args.distances)
a_groups, a_labels, a_data = read_long_csv(args.angles)

# Replicas — use intersection in the order they appear in the distance CSV
d_reps = list(d_data.keys())
a_reps = list(a_data.keys())
common_reps = [r for r in d_reps if r in a_reps]
missing_in_a = [r for r in d_reps if r not in a_reps]
missing_in_d = [r for r in a_reps if r not in d_reps]
if missing_in_a or missing_in_d:
    print(f"WARN: replica mismatch — distances={d_reps}, dihedrals={a_reps}")
    print(f"      using common: {common_reps}")


def unique(seq):
    seen = []
    for x in seq:
        if x not in seen:
            seen.append(x)
    return seen


d_unique_groups = unique(d_groups)
a_unique_groups = unique(a_groups)

# Rows: dihedral groups first, then distance groups
all_groups = a_unique_groups + d_unique_groups
n_rows = len(all_groups)
n_cols = len(common_reps)

a_group_to_row = {g: i + 1                       for i, g in enumerate(a_unique_groups)}
d_group_to_row = {g: i + 1 + len(a_unique_groups) for i, g in enumerate(d_unique_groups)}

print(f"Layout: {n_rows} rows ({len(a_unique_groups)} ang + {len(d_unique_groups)} dist) × {n_cols} cols")
print(f"        distance pairs: {len(d_labels)}, dihedral angles: {len(a_labels)}")


# ── Helpers ─────────────────────────────────────────────────────────
def smooth(y, window):
    """Running mean for linear data (distances)."""
    kernel = np.ones(window) / window
    return np.convolve(y, kernel, mode="same")


def break_wraps(t, y, threshold=180):
    """Insert NaN where |Δy| > threshold so Plotly breaks the line."""
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


# ── Plot ─────────────────────────────────────────────────────────────
fig = make_subplots(
    rows=n_rows, cols=n_cols,
    shared_xaxes="columns",  # all rows of a replica share its time axis
    shared_yaxes="rows",     # same group, same y across replicas
    vertical_spacing=0.0004, horizontal_spacing=0.04,
    row_titles=all_groups, column_titles=common_reps,
)

colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

# ─── Distance traces ───
for j, (label, group) in enumerate(zip(d_labels, d_groups)):
    row_idx = d_group_to_row[group]
    color = colors[j % len(colors)]
    for col_idx, rep in enumerate(common_reps, start=1):
        times, all_dists = d_data[rep]
        ref_d = all_dists[j, 0]
        # TODO: toggle initial value subtraction
        raw = all_dists[j] - ref_d
        # raw = all_dists[j]
        show_in_legend = (col_idx == 1)

        fig.add_trace(
            go.Scatter(
                x=times, y=raw, mode="lines", name=label,
                line=dict(width=0.5, color=color), opacity=0.3,
                showlegend=show_in_legend, legendgroup=label,
                hovertemplate=f"{label} ({rep})<br>t=%{{x:.2f}} ns<br>d=%{{y:.2f}} Å<extra></extra>",
            ),
            row=row_idx, col=col_idx,
        )
        if args.window > 1:
            fig.add_trace(
                go.Scatter(
                    x=times, y=smooth(raw, args.window), mode="lines",
                    name=f"{label} (avg {args.window}f)", showlegend=False,
                    line=dict(width=2, color=color), legendgroup=label,
                    hovertemplate=f"{label} ({rep}) smoothed<br>t=%{{x:.2f}} ns<br>d=%{{y:.2f}} Å<extra></extra>",
                ),
                row=row_idx, col=col_idx,
            )

# ─── Dihedral traces ───
# j resets to 0 inside this loop, so the dihedral palette starts fresh.
opa_raw = 0.3 if args.window > 1 else 1.0
for j, (label, group) in enumerate(zip(a_labels, a_groups)):
    row_idx = a_group_to_row[group]
    color = colors[j % len(colors)]
    for col_idx, rep in enumerate(common_reps, start=1):
        times, all_angles = a_data[rep]
        angles = all_angles[j]
        show_in_legend = (col_idx == 1)

        fig.add_trace(
            go.Scatter(
                x=times, y=angles, mode="markers", name=label,
                marker=dict(size=3, color=color), opacity=opa_raw,
                showlegend=show_in_legend, legendgroup=label,
                hovertemplate=f"{label} ({rep})<br>t=%{{x:.2f}} ns<br>θ=%{{y:.1f}}°<extra></extra>",
            ),
            row=row_idx, col=col_idx,
        )
        if args.window > 1:
            smoothed = circular_smooth(angles, args.window)
            t_sm, y_sm = break_wraps(times, smoothed, args.threshold)
            fig.add_trace(
                go.Scatter(
                    x=t_sm, y=y_sm, mode="lines",
                    name=f"{label} (avg {args.window}f)", showlegend=False,
                    line=dict(width=2, color=color), legendgroup=label,
                    hovertemplate=f"{label} ({rep}) smoothed<br>t=%{{x:.2f}} ns<br>θ=%{{y:.1f}}°<extra></extra>",
                ),
                row=row_idx, col=col_idx,
            )

# ── Axes ────────────────────────────────────────────────────────────
# Distance rows: pin every distance row to one global y-range so groups with
# small variations stay visible alongside groups with large variations
# (mirrors shared_yaxes="all" behaviour of the standalone distance plot).
all_dist_raw = []
for j in range(len(d_labels)):
    for rep in common_reps:
        _, all_dists = d_data[rep]
        all_dist_raw.append(all_dists[j] - all_dists[j, 0])
all_dist_raw = np.concatenate(all_dist_raw)
d_pad = (all_dist_raw.max() - all_dist_raw.min()) * 0.05
d_range = [all_dist_raw.min() - d_pad, all_dist_raw.max() + d_pad]

for g in d_unique_groups:
    fig.update_yaxes(title_text="Distance (Å)", range=d_range,
                     row=d_group_to_row[g], col=1)

# Dihedral rows: fixed range with rotamer ticks, "Angle (°)" label
for g in a_unique_groups:
    fig.update_yaxes(title_text="Angle (°)",
                     range=[-200, 200],
                     tickvals=[-180, -120, -60, 0, 60, 120, 180],
                     row=a_group_to_row[g], col=1)

# X-axis title on bottom row
for c in range(1, n_cols + 1):
    fig.update_xaxes(title_text="Time (ns)", row=n_rows, col=c)

fig.update_layout(
    height=max(900, 220 * n_rows),
    width=max(900, 450 * n_cols) + 280,  # extra width for the right-side legend
    template="plotly_white",
    legend=dict(orientation="v", yanchor="top", y=1.0, xanchor="left", x=1.02,
                font=dict(size=9)),
    title="Combined Distance + Dihedral Analysis",
)

out_html = f"{args.prefix}_combined_plots.html"
fig.write_html(out_html, include_plotlyjs="cdn")
print(f"Saved {out_html}")
fig.show()