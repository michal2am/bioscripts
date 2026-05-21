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
    # ("A:Phe200-CA ↔ B:Phe46-CA", "A:Lys197-NZ ↔ A:Glu165-OE1", "xxx"),
    # ("C:Phe200-CA ↔ D:Phe46-CA", "C:Lys197-NZ ↔ C:Glu165-OE1", "yyy"),
    # ("A:Phe200-CA ↔ B:Phe46-CA", "A:Lys196-NZ ↔ A:Glu153-OE1", "xxx"),
    # ("C:Phe200-CA ↔ D:Phe46-CA", "C:Lys196-NZ ↔ C:Glu153-OE1", "yyy"),

    ("A:X195-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    ("C:X195-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),

    ("A:Lys196-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    ("C:Lys196-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),

    ("A:X197-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    ("C:X197-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),

    ("A:X198-CA ↔ B:Phe46-CA", "A:Phe200-CA ↔ B:Phe46-CA", "BS_1"),
    ("C:X198-CA ↔ D:Phe46-CA", "C:Phe200-CA ↔ D:Phe46-CA", "BS_2"),
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
    # shared_xaxes="rows",  # same scatter pair → same x range across replicas
    # shared_yaxes="rows",  # same scatter pair → same y range across replicas
    shared_xaxes="all",
    shared_yaxes="all",
    vertical_spacing=0.04, horizontal_spacing=0.04,
    row_titles=[t for _, _, _, _, t in resolved],
    column_titles=unique_replicas,
)

# Global time range so the viridis mapping is consistent across all panels —
# late frames in a long replica and late frames in a short replica still get
# the same colour at the same absolute time.
all_times = np.concatenate([replica_data[rep][0] for rep in unique_replicas])
t_min, t_max = all_times.min(), all_times.max()

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
            go.Scatter(
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

            fig.add_trace(
                go.Scatter(
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