#!/usr/bin/env python3
"""Plot per-residue dihedral angles from a CSV produced by dihedrals_calc.py.

Raw data is shown as scatter markers (no transformations), with a circular-mean
smoothing line on top. The smoothed line uses NaN-break at ±180° wraps to avoid
drawing fake ~360° jumps.
"""

import argparse
import csv
import numpy as np
import plotly.graph_objects as go
from astropy.units.format import ogip_parsetab
from plotly.subplots import make_subplots

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="Plot dihedral angles from CSV")
p.add_argument("-i", "--input", required=True, help="Input CSV (from dihedrals_calc.py)")
p.add_argument("-w", "--window", type=int, default=50, help="Smoothing window in frames (default: 50, 0=off)")
p.add_argument("-t", "--threshold", type=float, default=180,
               help="Wrap-detection threshold in deg for the smoothed line (default: 180)")
p.add_argument("-o", "--prefix", default="dihedrals", help="Output file prefix")
args = p.parse_args()

# ── Read CSV ─────────────────────────────────────────────────────────
with open(args.input) as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = list(reader)

data = np.array(rows, dtype=float)
times = data[:, 0]
all_angles = data[:, 1:].T  # (n_series, n_frames), in degrees

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

print(f"Loaded {len(times)} frames, {len(labels)} angles, {len(unique_groups)} groups: {unique_groups}")


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


# ── Plot ─────────────────────────────────────────────────────────────
n_rows = len(unique_groups)
fig = make_subplots(rows=n_rows, cols=1, shared_xaxes=True, vertical_spacing=0.06,
                    subplot_titles=tuple(unique_groups))

colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

for j, (label, group) in enumerate(zip(labels, groups)):
    row = group_to_row[group]
    angles = all_angles[j]
    color = colors[j % len(colors)]

    # Raw scatter (markers only — no line so no wrap artefacts)
    if args.window > 1:
        opa=0.3
    else :
        opa=1

    fig.add_trace(
        go.Scatter(
            x=times, y=angles, mode="markers", name=label,
            marker=dict(size=3, color=color), opacity=opa,
            hovertemplate=f"{label}<br>t=%{{x:.2f}} ns<br>θ=%{{y:.1f}}°<extra></extra>",
            legendgroup=label,
        ),
        row=row, col=1,
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
                hovertemplate=f"{label} smoothed<br>t=%{{x:.2f}} ns<br>θ=%{{y:.1f}}°<extra></extra>",
                legendgroup=label,
            ),
            row=row, col=1,
        )

for r in range(1, n_rows + 1):
    fig.update_yaxes(title_text="Angle (°)", row=r, col=1,
                     range=[-200, 200],
                     tickvals=[-180, -120, -60, 0, 60, 120, 180])
fig.update_xaxes(title_text="Time (ns)", row=n_rows, col=1)

fig.update_layout(
    height=900, template="plotly_white",
    legend=dict(orientation="h", yanchor="bottom", y=1.03, xanchor="center", x=0.5,
                font=dict(size=9)),
    title="Dihedral Angles",
)

out_html = f"{args.prefix}_dihedrals_plots.html"
fig.write_html(out_html, include_plotlyjs="cdn")
print(f"Saved {out_html}")
fig.show()