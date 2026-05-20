#!/usr/bin/env python3
"""
Plot TICA results produced by tica_calc.py.

Reads <prefix>_tics.csv, <prefix>_eigenvalues.csv, <prefix>_loadings.csv
and writes an interactive HTML figure.
"""

import argparse
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


p = argparse.ArgumentParser(description="Plot TICA results")
p.add_argument("-i", "--input-prefix", required=True, help="Input prefix (from tica_calc.py)")
p.add_argument("-o", "--output", default=None, help="Output HTML file (default: <prefix>_plot.html)")
p.add_argument("--n-tics-shown", type=int, default=2, help="Number of TICs in time series (default: 2)")
p.add_argument("--fel-bins", type=int, default=60, help="Bins for free energy landscape (default: 60)")
p.add_argument("--fel-max", type=float, default=6.0, help="Cap for −ln P colorscale in kT (default: 6.0)")
args = p.parse_args()


# Load
tics = pd.read_csv(f"{args.input_prefix}_tics.csv")
eig = pd.read_csv(f"{args.input_prefix}_eigenvalues.csv")
loadings = pd.read_csv(f"{args.input_prefix}_loadings.csv")

# Discover sites
sites = []
for col in tics.columns:
    if "|" in col:
        s = col.split("|")[0]
        if s not in sites:
            sites.append(s)

n_sites = len(sites)
n_show = args.n_tics_shown
n_rows = n_show + 2
fel_row = n_show + 1
load_row = n_show + 2

# Subplot titles (row-major)
titles = []
for i in range(n_show):
    for s in sites:
        titles.append(f"{s}: TIC{i+1} vs time")
for s in sites:
    titles.append(f"{s}: free energy landscape")
for s in sites:
    titles.append(f"{s}: feature ↔ TIC correlations")

fig = make_subplots(
    rows=n_rows, cols=n_sites,
    subplot_titles=titles,
    vertical_spacing=0.08,
    horizontal_spacing=0.12,
)

# --- TIC time series ---
tic_colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b"]
for tic_i in range(n_show):
    for col_i, site in enumerate(sites):
        col_name = f"{site}|TIC{tic_i+1}"
        fig.add_trace(
            go.Scattergl(
                x=tics["time_ns"], y=tics[col_name],
                mode="lines",
                line=dict(width=1, color=tic_colors[tic_i % len(tic_colors)]),
                showlegend=False, name=col_name,
            ),
            row=tic_i+1, col=col_i+1
        )

# --- Free energy landscape (shared scale across sites) ---
all_x = np.concatenate([tics[f"{s}|TIC1"].values for s in sites])
all_y = np.concatenate([tics[f"{s}|TIC2"].values for s in sites])
xe = np.linspace(all_x.min(), all_x.max(), args.fel_bins + 1)
ye = np.linspace(all_y.min(), all_y.max(), args.fel_bins + 1)
xc = 0.5 * (xe[:-1] + xe[1:])
yc = 0.5 * (ye[:-1] + ye[1:])

for col_i, site in enumerate(sites):
    x = tics[f"{site}|TIC1"].values
    y = tics[f"{site}|TIC2"].values
    H, _, _ = np.histogram2d(x, y, bins=[xe, ye], density=True)
    H = np.where(H > 0, H, np.nan)
    F = -np.log(H)
    F = F - np.nanmin(F)
    F = np.clip(F, 0, args.fel_max)
    fig.add_trace(
        go.Heatmap(x=xc, y=yc, z=F.T, coloraxis="coloraxis"),
        row=fel_row, col=col_i+1
    )

# --- Loadings heatmap (shared [-1, 1]) ---
features = loadings["feature"].tolist()
for col_i, site in enumerate(sites):
    site_cols = [c for c in loadings.columns if c.startswith(f"{site}|TIC")]
    L = loadings[site_cols].values
    tic_labels = [c.split("|")[1] for c in site_cols]
    fig.add_trace(
        go.Heatmap(
            x=tic_labels, y=features, z=L,
            text=np.round(L, 2), texttemplate="%{text}",
            textfont=dict(size=10),
            coloraxis="coloraxis2",
        ),
        row=load_row, col=col_i+1
    )

# --- Axis labels ---
for tic_i in range(n_show):
    for col_i in range(n_sites):
        fig.update_yaxes(title_text=f"TIC{tic_i+1}", row=tic_i+1, col=col_i+1)
        if tic_i == n_show - 1:
            fig.update_xaxes(title_text="time [ns]", row=tic_i+1, col=col_i+1)

for col_i in range(n_sites):
    fig.update_xaxes(title_text="TIC1", row=fel_row, col=col_i+1)
    fig.update_yaxes(title_text="TIC2", row=fel_row, col=col_i+1)
    fig.update_xaxes(title_text="TIC", row=load_row, col=col_i+1)
    fig.update_yaxes(title_text="H-bond", row=load_row, col=col_i+1, tickmode="linear", dtick=1,)

# --- Eigenvalues / timescales in subtitle ---
n_tic_avail = len(eig)
subtitle_lines = []
for s in sites:
    eigs = eig[f"{s}|eigenvalue"].values
    ts = eig[f"{s}|timescale_ps"].values
    parts = [f"TIC{i+1}: λ={eigs[i]:.3f}, τ={ts[i]:.0f} ps"
             for i in range(min(4, n_tic_avail))]
    subtitle_lines.append(f"<b>{s}</b>  " + " | ".join(parts))
subtitle = "<br>".join(subtitle_lines)

# --- Colorbar positions tied to row centers ---
def row_y_center(row):
    return 1.0 - (row - 0.5) / n_rows

fig.update_layout(
    title=dict(text=f"TICA: B9-B10 H-bonds<br><sub>{subtitle}</sub>", x=0.5),
    height=320 * n_rows,
    width=560 * n_sites + 220,
    template="plotly_white",
    coloraxis=dict(
        colorscale="Viridis",
        cmin=0, cmax=args.fel_max,
        colorbar=dict(title="−ln P", x=1.02,
                      y=row_y_center(fel_row),
                      len=0.85/n_rows, yanchor="middle"),
    ),
    coloraxis2=dict(
        colorscale="RdBu",
        cmin=-1, cmax=1,
        colorbar=dict(title="r", x=1.02,
                      y=row_y_center(load_row),
                      len=0.85/n_rows, yanchor="middle"),
    ),
)

output = args.output or f"{args.input_prefix}_tics_plot.html"
fig.write_html(output)
print(f"Wrote {output}")