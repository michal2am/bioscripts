#!/usr/bin/env python3
"""Plot atom-pair distances from a CSV produced by distances_calc.py."""

import argparse
import csv
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="Plot atom-pair distances from CSV")
p.add_argument("-i", "--input", required=True, help="Input CSV (from distances_calc.py)")
p.add_argument("-w", "--window", type=int, default=50, help="Smoothing window in frames (default: 50, 0=off)")
p.add_argument("-o", "--prefix", default="distances", help="Output file prefix")
args = p.parse_args()

# ── Read CSV ─────────────────────────────────────────────────────────
with open(args.input) as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = list(reader)

data = np.array(rows, dtype=float)
times = data[:, 0]
all_dists = data[:, 1:].T  # (n_pairs, n_frames)

# Each column header is "group|label"
columns = header[1:]
groups = [c.split("|", 1)[0] for c in columns]
labels = [c.split("|", 1)[1] for c in columns]

# Unique groups in order of first appearance → one subplot row each
unique_groups = []
for g in groups:
    if g not in unique_groups:
        unique_groups.append(g)
group_to_row = {g: i + 1 for i, g in enumerate(unique_groups)}

print(f"Loaded {len(times)} frames, {len(labels)} pairs, {len(unique_groups)} groups: {unique_groups}")


# ── Interactive Plotly figure ────────────────────────────────────────
def smooth(y, window=50):
    """Running average with edge handling."""
    kernel = np.ones(window) / window
    return np.convolve(y, kernel, mode="same")


n_rows = len(unique_groups)
fig = make_subplots(rows=n_rows, cols=1, shared_xaxes=True, vertical_spacing=0.06,
                    subplot_titles=tuple(unique_groups))

colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

for j, (label, group) in enumerate(zip(labels, groups)):
    row = group_to_row[group]
    ref_d = all_dists[j, 0]
    #TODO: toggle initial value substraction
    raw = all_dists[j] - ref_d
    #raw = all_dists[j]
    color = colors[j % len(colors)]

    # Raw trace (thin, semi-transparent)
    fig.add_trace(
        go.Scatter(
            x=times, y=raw, mode="lines", name=label,
            line=dict(width=0.5, color=color), opacity=0.3,
            hovertemplate=f"{label}<br>t=%{{x:.2f}} ns<br>d=%{{y:.2f}} Å<extra></extra>",
            legendgroup=label,
        ),
        row=row, col=1,
    )
    # Smoothed trace (bold, on top)
    if args.window > 1:
        fig.add_trace(
            go.Scatter(
                x=times, y=smooth(raw, args.window), mode="lines",
                name=f"{label} (avg {args.window}f)", showlegend=False,
                line=dict(width=2, color=color),
                hovertemplate=f"{label} smoothed<br>t=%{{x:.2f}} ns<br>d=%{{y:.2f}} Å<extra></extra>",
                legendgroup=label,
            ),
            row=row, col=1,
        )
    # add initial-frame reference as dashed line (toggles with main trace)
    '''
    color = colors[j % len(colors)]
    fig.add_trace(
        go.Scatter(
            x=[times[0], times[-1]], y=[ref_d, ref_d], mode="lines+text",
            line=dict(width=1, dash="dot", color=color),
            opacity=0.5, showlegend=False, legendgroup=label,
            text=[f"t₀={ref_d:.2f} Å", ""], textposition="top left",
            textfont=dict(size=10, color=color),
            hoverinfo="skip",
        ),
        row=row, col=1,

    )
    '''

for r in range(1, n_rows + 1):
    fig.update_yaxes(title_text="Distance (Å)", row=r, col=1)
fig.update_xaxes(title_text="Time (ns)", row=n_rows, col=1)

fig.update_layout(
    height=900, template="plotly_white",
    legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="center", x=0.5,
                font=dict(size=9)),
    title="Atom-Pair Distances",
)

out_html = f"{args.prefix}_distances_plots.html"
fig.write_html(out_html, include_plotlyjs="cdn")
print(f"Saved {out_html}")
fig.show()