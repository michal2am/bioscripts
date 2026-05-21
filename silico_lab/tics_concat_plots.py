#!/usr/bin/env python3
"""
Plot TICA results from tica_calc.py outputs (multi-replica aware).

Layout:
  - 1 row (optional): auxiliary interface distances per site (from <prefix>_aux.csv)
  - n_tics-shown rows: TIC time series, one trace per replica (color = replica)
  - n_replicas rows:   TIC1 vs TIC2 free energy landscape, one row per replica
  - 1 row:             feature ↔ TIC correlation heatmaps
  Columns = binding sites.
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
p.add_argument("--loadings-height", type=float, default=1.0,
               help="Loadings row height as multiplier of other rows (default: 1.0)")
args = p.parse_args()


# Load
tics = pd.read_csv(f"{args.input_prefix}_tics.csv")
eig = pd.read_csv(f"{args.input_prefix}_eigenvalues.csv")
loadings = pd.read_csv(f"{args.input_prefix}_loadings.csv")

# Optional auxiliary distances (interface contacts tracked alongside TICA)
try:
    aux = pd.read_csv(f"{args.input_prefix}_aux.csv")
    has_aux = True
except FileNotFoundError:
    has_aux = False

# Discover sites
sites = []
for col in tics.columns:
    if "|" in col:
        s = col.split("|")[0]
        if s not in sites:
            sites.append(s)
n_sites = len(sites)

# Discover replicas (backwards-compatible with single-replica CSVs)
if "replica" in tics.columns:
    replicas = sorted(tics["replica"].unique())
else:
    tics["replica"] = 1
    replicas = [1]
n_rep = len(replicas)

n_show = args.n_tics_shown
n_aux = 1 if has_aux else 0
aux_row = 1
ts_first_row = n_aux + 1
fel_first_row = n_aux + n_show + 1
fel_last_row = n_aux + n_show + n_rep
load_row = n_aux + n_show + n_rep + 1
n_rows = load_row

# Subplot titles (row-major)
titles = []
if has_aux:
    for s in sites:
        aux_col = next(c for c in aux.columns if c.startswith(f"{s}|"))
        titles.append(f"{s}: {aux_col.split('|')[1]}")
for i in range(n_show):
    for s in sites:
        titles.append(f"{s}: TIC{i+1} vs time")
for rep in replicas:
    for s in sites:
        titles.append(f"{s}: FEL (replica {rep})")
for s in sites:
    titles.append(f"{s}: feature ↔ TIC correlations")

row_heights = [1.0] * (n_rows - 1) + [args.loadings_height]

fig = make_subplots(
    rows=n_rows, cols=n_sites,
    subplot_titles=titles,
    vertical_spacing=0.05,
    horizontal_spacing=0.12,
    row_heights=row_heights,
)

replica_colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b"]

# --- Auxiliary distance time series (top row, if present) ---
if has_aux:
    for rep_i, rep in enumerate(replicas):
        rep_data = aux[aux["replica"] == rep]
        color = replica_colors[rep_i % len(replica_colors)]
        for col_i, site in enumerate(sites):
            aux_col = next(c for c in aux.columns if c.startswith(f"{site}|"))
            fig.add_trace(
                go.Scattergl(
                    x=rep_data["time_ns"], y=rep_data[aux_col],
                    mode="lines",
                    line=dict(width=1, color=color),
                    name=f"replica {rep}",
                    legendgroup=f"rep{rep}",
                    showlegend=(col_i == 0),
                ),
                row=aux_row, col=col_i+1
            )

# --- TIC time series, one trace per replica per panel ---
for rep_i, rep in enumerate(replicas):
    rep_data = tics[tics["replica"] == rep]
    color = replica_colors[rep_i % len(replica_colors)]
    for tic_i in range(n_show):
        for col_i, site in enumerate(sites):
            col_name = f"{site}|TIC{tic_i+1}"
            fig.add_trace(
                go.Scattergl(
                    x=rep_data["time_ns"], y=rep_data[col_name],
                    mode="lines",
                    line=dict(width=1, color=color),
                    name=f"replica {rep}",
                    legendgroup=f"rep{rep}",
                    showlegend=(not has_aux and tic_i == 0 and col_i == 0),
                ),
                row=ts_first_row + tic_i, col=col_i+1
            )

# --- Free energy landscapes, one row per replica (shared bins + color scale) ---
all_x = np.concatenate([tics[f"{s}|TIC1"].values for s in sites])
all_y = np.concatenate([tics[f"{s}|TIC2"].values for s in sites])
xe = np.linspace(all_x.min(), all_x.max(), args.fel_bins + 1)
ye = np.linspace(all_y.min(), all_y.max(), args.fel_bins + 1)
xc = 0.5 * (xe[:-1] + xe[1:])
yc = 0.5 * (ye[:-1] + ye[1:])

for rep_i, rep in enumerate(replicas):
    rep_data = tics[tics["replica"] == rep]
    for col_i, site in enumerate(sites):
        x = rep_data[f"{site}|TIC1"].values
        y = rep_data[f"{site}|TIC2"].values
        H, _, _ = np.histogram2d(x, y, bins=[xe, ye], density=True)
        H = np.where(H > 0, H, np.nan)
        F = -np.log(H)
        F = F - np.nanmin(F)
        F = np.clip(F, 0, args.fel_max)
        fig.add_trace(
            go.Heatmap(x=xc, y=yc, z=F.T, coloraxis="coloraxis"),
            row=fel_first_row + rep_i, col=col_i+1
        )

# --- Loadings heatmap (shared [-1, 1]) ---
features = loadings["feature"].tolist()
for col_i, site in enumerate(sites):
    site_cols = [c for c in loadings.columns if c.startswith(f"{site}|TIC")][:10]
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
if has_aux:
    for col_i in range(n_sites):
        fig.update_yaxes(title_text="d [Å]", row=aux_row, col=col_i+1)

for tic_i in range(n_show):
    for col_i in range(n_sites):
        fig.update_yaxes(title_text=f"TIC{tic_i+1}", row=ts_first_row + tic_i, col=col_i+1)
        if tic_i == n_show - 1:
            fig.update_xaxes(title_text="time [ns]", row=ts_first_row + tic_i, col=col_i+1)

for rep_i in range(n_rep):
    for col_i in range(n_sites):
        fig.update_xaxes(title_text="TIC1", row=fel_first_row + rep_i, col=col_i+1)
        fig.update_yaxes(title_text="TIC2", row=fel_first_row + rep_i, col=col_i+1)

for col_i in range(n_sites):
    fig.update_xaxes(title_text="TIC", row=load_row, col=col_i+1)
    fig.update_yaxes(title_text="H-bond", tickmode="linear", dtick=1, row=load_row, col=col_i+1)

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

# --- Colorbar positions (account for non-uniform row heights) ---
total_h = sum(row_heights)
cum = [0.0]
for h in row_heights:
    cum.append(cum[-1] + h / total_h)

fel_top = cum[fel_first_row - 1]
fel_bot = cum[fel_last_row]
fel_y = 1.0 - 0.5 * (fel_top + fel_bot)
fel_len = 0.85 * (fel_bot - fel_top)

load_top = cum[load_row - 1]
load_bot = cum[load_row]
load_y = 1.0 - 0.5 * (load_top + load_bot)
load_len = 0.85 * (load_bot - load_top)

fig.update_layout(
    title=dict(text=f"<sub>{subtitle}</sub>", x=0.5),
    height=int(300 * total_h),
    width=500 * n_sites + 220,
    template="plotly_white",
    legend=dict(orientation="h", x=0.5, xanchor="center", y=1.0, yanchor="bottom"),
    coloraxis=dict(
        colorscale="Viridis",
        cmin=0, cmax=args.fel_max,
        colorbar=dict(title="−ln P", x=1.02, y=fel_y, len=fel_len, yanchor="middle"),
    ),
    coloraxis2=dict(
        colorscale="RdBu",
        cmin=-1, cmax=1,
        colorbar=dict(title="r", x=1.02, y=load_y, len=load_len, yanchor="middle"),
    ),
)

output = args.output or f"{args.input_prefix}_tics_plot.html"
fig.write_html(output)
print(f"Wrote {output}")