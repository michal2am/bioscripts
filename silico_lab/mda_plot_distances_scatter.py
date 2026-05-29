#!/usr/bin/env python3
"""Plot x-y scatters of selected distance pairs from a CSV produced by
distances_calc.py.

Layout: rows = scatter pairs (defined in SCATTER_PAIRS), columns = replicas.
Each frame is one point, colored by simulation time so the trajectory through
state space is visible.
"""

import argparse
import csv
from collections import defaultdict
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ── Scatter pair definitions ────────────────────────────────────────
# Each entry: (x_label, y_label, panel_title)
# Labels must match the "label" portion of the CSV column headers.
SCATTER_PAIRS = [

    #("A:Phe200-CA ↔ B:Phe46-CA", "A:Tyr205-OH ↔ A:Thr202-OG1", "BS_1"),
    #("A:Phe200-CA ↔ B:Phe46-CA", "A:Asp162-CG ↔ A:Thr202-OG1", "BS_1"),

    #("C:Phe200-CA ↔ D:Phe46-CA", "C:Tyr205-OH ↔ C:Thr202-OG1", "BS_2"),
    #("C:Phe200-CA ↔ D:Phe46-CA", "C:Asp162-CG ↔ C:Thr202-OG1", "BS_2"),


    # ("A:Phe200-CA ↔ B:Phe46-CA", "A:Lys197-NZ ↔ A:Glu165-OE1", "xxx"),
    # ("C:Phe200-CA ↔ D:Phe46-CA", "C:Lys197-NZ ↔ C:Glu165-OE1", "yyy"),
    # ("A:Phe200-CA ↔ B:Phe46-CA", "A:Lys196-NZ ↔ A:Glu153-OE1", "xxx"),
    # ("C:Phe200-CA ↔ D:Phe46-CA", "C:Lys196-NZ ↔ C:Glu153-OE1", "yyy"),
    #
    # ('A:glyc2-O ↔ A:Lys196-NZ', "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    # ('C:glyc2-O ↔ C:Lys196-NZ', "C:Phe200-CA ↔ C:Phe46-CA", "BS_2"),

    # ('A:glyc2-O ↔ A:Lys196-NZ', "A:Lys196-NZ ↔ A:Glu153-OE1", "BS_1"),
    # ('C:glyc2-O ↔ C:Lys196-NZ', "C:Lys196-NZ ↔ C:Glu153-OE1", "BS_2"),


    # loop C vs F46

    # ("A:X194-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    # ("C:X194-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),

    ("A:X195-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    ("C:X195-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),

    # ("A:Lys196-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    # ("C:Lys196-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),
    #
    # ("A:X197-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    # ("C:X197-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),
    #
    # ("A:X198-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    # ("C:X198-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),

    ("A:X199-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    ("C:X199-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),

    ("A:X201-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    ("C:X201-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),

    # ("A:X202-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    # ("C:X202-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),
    #
    # ("A:X203-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    # ("C:X203-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),
    #
    # ("A:X204-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    # ("C:X204-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),

    ("A:X205-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    ("C:X205-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),

]

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="Plot x-y distance scatters from CSV")
p.add_argument("-i", "--input", required=True, help="Input CSV (from distances_calc.py)")
p.add_argument("-o", "--prefix", default="distances", help="Output file prefix")
p.add_argument("--marker-size", type=float, default=3.0, help="Marker size (default: 3.0)")
p.add_argument("--alpha", type=float, default=0.5, help="Marker opacity (default: 0.5)")
args = p.parse_args()

# ── Read CSV ─────────────────────────────────────────────────────────
# Header: time_ns, replica, "group|label", "group|label", ...
with open(args.input) as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = list(reader)

columns = header[2:]
labels_all = [c.split("|", 1)[1] for c in columns]
label_to_idx = {l: i for i, l in enumerate(labels_all)}

# Group rows by replica, preserving insertion order
rep_rows = defaultdict(list)
for row in rows:
    rep_rows[row[1]].append(row)

unique_replicas = list(rep_rows.keys())

# Per-replica numpy arrays
replica_data = {}
for rep, rrows in rep_rows.items():
    times = np.array([float(r[0]) for r in rrows])
    dists = np.array([[float(x) for x in r[2:]] for r in rrows]).T
    replica_data[rep] = (times, dists)

# Resolve each scatter pair's column indices against the CSV labels
resolved = []
for x_lbl, y_lbl, title in SCATTER_PAIRS:
    if x_lbl not in label_to_idx or y_lbl not in label_to_idx:
        missing = x_lbl if x_lbl not in label_to_idx else y_lbl
        print(f"  SKIP  '{title}' — label not in CSV: {missing}")
        continue
    resolved.append((label_to_idx[x_lbl], label_to_idx[y_lbl], x_lbl, y_lbl, title))
    print(f"  OK    {title}")

n_rows = len(resolved)
n_cols = len(unique_replicas)

print(f"Layout: {n_rows} scatter pairs × {n_cols} replicas: {unique_replicas}")

# ── Plot ─────────────────────────────────────────────────────────────
fig = make_subplots(
    rows=n_rows, cols=n_cols,
    shared_xaxes="rows",  # same scatter pair → same x range across replicas
    shared_yaxes="rows",  # same scatter pair → same y range across replicas
    vertical_spacing=0.008, horizontal_spacing=0.04,
    row_titles=[t for _, _, _, _, t in resolved],
    column_titles=unique_replicas,
)

# Global time range so the viridis mapping is consistent across all panels —
# late frames in a long replica and late frames in a short replica still get
# the same colour at the same absolute time.
all_times = np.concatenate([replica_data[rep][0] for rep in unique_replicas])
t_min, t_max = all_times.min(), all_times.max()

# Collect R² per (replica, title, x_lbl) for the summary figure built after
r2_data = {}

for r_idx, (ix, iy, x_lbl, y_lbl, title) in enumerate(resolved):
    row = r_idx + 1
    for c_idx, rep in enumerate(unique_replicas):
        col = c_idx + 1
        times, all_dists = replica_data[rep]
        # TODO: toggle initial value subtraction
        x_vals = all_dists[ix] - all_dists[ix, 0]
        y_vals = all_dists[iy] - all_dists[iy, 0]
        # x_vals = all_dists[ix]
        # y_vals = all_dists[iy]

        # Only one trace draws the colorbar (top-right corner)
        show_colorbar = (r_idx == 0 and c_idx == n_cols - 1)

        fig.add_trace(
            go.Scattergl(
                x=x_vals, y=y_vals, mode="markers",
                marker=dict(
                    size=args.marker_size,
                    color=times,
                    colorscale="Viridis",
                    cmin=t_min, cmax=t_max,
                    showscale=show_colorbar,
                    colorbar=dict(title="Time (ns)") if show_colorbar else None,
                ),
                opacity=args.alpha,
                showlegend=False,
                hovertemplate=(
                    f"{title} ({rep})<br>"
                    f"x=%{{x:.2f}} Å<br>y=%{{y:.2f}} Å<extra></extra>"
                ),
            ),
            row=row, col=col,
        )

        # Linear regression for this panel
        if len(x_vals) >= 2 and np.std(x_vals) > 0:
            slope, intercept = np.polyfit(x_vals, y_vals, 1)
            x_line = np.array([x_vals.min(), x_vals.max()])
            y_line = slope * x_line + intercept
            y_pred = slope * x_vals + intercept
            ss_res = np.sum((y_vals - y_pred) ** 2)
            ss_tot = np.sum((y_vals - y_vals.mean()) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
            r2_data[(rep, title, x_lbl)] = r2

            fig.add_trace(
                go.Scattergl(
                    x=x_line, y=y_line, mode="lines",
                    line=dict(width=2, color="crimson"),
                    showlegend=False,
                    hovertemplate=(
                        f"{title} ({rep}) regression<br>"
                        f"y = {slope:.3f}·x + {intercept:.3f}<br>"
                        f"R² = {r2:.3f}<extra></extra>"
                    ),
                ),
                row=row, col=col,
            )
            fig.add_annotation(
                text=f"R² = {r2:.2f}, slope = {slope:.2f}",
                xref="x domain", yref="y domain",
                x=0.04, y=0.96, xanchor="left", yanchor="top",
                showarrow=False,
                font=dict(size=10, color="crimson"),
                bgcolor="rgba(255,255,255,0.75)",
                bordercolor="crimson", borderwidth=0.5, borderpad=2,
                row=row, col=col,
            )

# Axis titles: y on leftmost column of each row; x on every subplot of each row
# (each row has its own distance pair on x, so it gets its own x label).
for r_idx, (ix, iy, x_lbl, y_lbl, title) in enumerate(resolved):
    row = r_idx + 1
    fig.update_yaxes(title_text=f"{y_lbl} (Å)", row=row, col=1)
    for c in range(1, n_cols + 1):
        fig.update_xaxes(title_text=f"{x_lbl} (Å)", row=row, col=c)

fig.update_layout(
    height=max(900, 400 * n_rows),  # scatter panels look better tall-ish
    width=max(900, 400 * n_cols) + 160,  # extra width for colorbar
    template="plotly_white",
    title="Atom-Pair Distance Scatters",
)

out_html = f"{args.prefix}_distances_scatter.html"
fig.write_html(out_html, include_plotlyjs="cdn")
print(f"Saved {out_html}")
fig.show()

# ── R² summary figure ───────────────────────────────────────────────
# Rows = binding-site titles, columns = replicas. Each panel shows R² of
# the linear regression for every x-label that belonged to that title.
unique_titles = []
title_to_xlbls = {}
for ix, iy, x_lbl, y_lbl, title in resolved:
    if title not in unique_titles:
        unique_titles.append(title)
        title_to_xlbls[title] = []
    if x_lbl not in title_to_xlbls[title]:
        title_to_xlbls[title].append(x_lbl)

n_r2_rows = len(unique_titles)
n_r2_cols = len(unique_replicas)

fig_r2 = make_subplots(
    rows=n_r2_rows, cols=n_r2_cols,
    shared_xaxes="rows",   # same binding site → same x labels across replicas
    shared_yaxes="all",    # R² is bounded [0, 1]
    vertical_spacing=0.16, horizontal_spacing=0.04,
    row_titles=unique_titles,
    column_titles=unique_replicas,
)

for r_idx, title in enumerate(unique_titles):
    row = r_idx + 1
    x_lbls = title_to_xlbls[title]
    for c_idx, rep in enumerate(unique_replicas):
        col = c_idx + 1
        r2_values = [r2_data.get((rep, title, x_lbl), float("nan")) for x_lbl in x_lbls]
        fig_r2.add_trace(
            go.Bar(
                x=x_lbls, y=r2_values,
                text=[f"{v:.2f}" if not np.isnan(v) else "" for v in r2_values],
                textposition="outside",
                marker=dict(color="steelblue"),
                showlegend=False,
                hovertemplate="%{x}<br>R² = %{y:.3f}<extra></extra>",
            ),
            row=row, col=col,
        )

# y-axis: fixed 0–1 range, "R²" title on leftmost column
for r in range(1, n_r2_rows + 1):
    fig_r2.update_yaxes(title_text="R²", range=[0, 1.05], row=r, col=1)

# x-axis: tilt labels (they can be long) on every panel
for r in range(1, n_r2_rows + 1):
    for c in range(1, n_r2_cols + 1):
        fig_r2.update_xaxes(tickangle=-30, row=r, col=c)

fig_r2.update_layout(
    height=max(600, 350 * n_r2_rows),
    width=max(900, 400 * n_r2_cols) + 160,
    template="plotly_white",
    title="Linear regression R² by binding site",
)

out_html_r2 = f"{args.prefix}_distances_scatter_r2.html"
fig_r2.write_html(out_html_r2, include_plotlyjs="cdn")
print(f"Saved {out_html_r2}")
fig_r2.show()