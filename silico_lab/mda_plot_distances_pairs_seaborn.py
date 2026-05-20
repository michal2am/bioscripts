#!/usr/bin/env python3
"""Publication-quality 2D scatter plots between pairs of distance series.

Reads the CSV from distances_calc.py and produces a side-by-side scatter
figure for the SCATTER_PAIRS listed below. Each point is one frame; points
are coloured by simulation time so the trajectory through the (x, y) state
space is visible.
"""

import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ── Scatter pair definitions ────────────────────────────────────────
# Each entry: (x_label, y_label, panel_title)
# Labels must match the "label" portion of CSV column headers.
SCATTER_PAIRS = [
    ("A:Tyr205-OH ↔ A:Thr202-OG1", "A:200-CG ↔ B:46-CG", "Chain A"),
    ("C:Tyr205-OH ↔ C:Thr202-OG1", "C:200-CG ↔ D:46-CG", "Chain C"),
]

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="2D distance scatter plots from CSV")
p.add_argument("-i", "--input", required=True, help="Input CSV (from distances_calc.py)")
p.add_argument("-o", "--prefix", default="distances", help="Output file prefix")
p.add_argument("--format", default="pdf", choices=["pdf", "png", "svg"], help="Output format (default: pdf)")
p.add_argument("--dpi", type=int, default=300, help="DPI for raster output (default: 300)")
p.add_argument("--width", type=float, default=7.0, help="Figure width in inches (default: 7.0)")
p.add_argument("--height", type=float, default=3.5, help="Figure height in inches (default: 3.5)")
p.add_argument("--marker-size", type=float, default=4.0, help="Marker size (default: 4.0)")
p.add_argument("--alpha", type=float, default=0.5, help="Marker alpha (default: 0.5)")
args = p.parse_args()

# ── Read CSV ─────────────────────────────────────────────────────────
with open(args.input) as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = list(reader)

data = np.array(rows, dtype=float)
times = data[:, 0]
all_dists = data[:, 1:].T  # (n_pairs, n_frames)

# Each column header is "group|label" — strip the group prefix
columns = header[1:]
labels = [c.split("|", 1)[1] for c in columns]
label_to_idx = {l: i for i, l in enumerate(labels)}

print(f"Loaded {len(times)} frames, {len(labels)} pairs")

# Resolve each scatter pair against the CSV columns
resolved = []
for x_lbl, y_lbl, title in SCATTER_PAIRS:
    if x_lbl not in label_to_idx or y_lbl not in label_to_idx:
        missing = x_lbl if x_lbl not in label_to_idx else y_lbl
        print(f"  SKIP  '{title}' — label not in CSV: {missing}")
        continue
    resolved.append((label_to_idx[x_lbl], label_to_idx[y_lbl], x_lbl, y_lbl, title))
    print(f"  OK    {title}")

# ── Styling ──────────────────────────────────────────────────────────
sns.set_theme(context="paper", style="ticks", font_scale=1.0)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# ── Figure ──────────────────────────────────────────────────────────
n_panels = len(resolved)
fig, axes = plt.subplots(1, n_panels,
                         figsize=(args.width, args.height),
                         constrained_layout=True)
axes = np.atleast_1d(axes)

# Shared axis ranges so panels are directly comparable
x_all = np.concatenate([all_dists[ix] for ix, _, _, _, _ in resolved])
y_all = np.concatenate([all_dists[iy] for _, iy, _, _, _ in resolved])
x_pad = (x_all.max() - x_all.min()) * 0.03
y_pad = (y_all.max() - y_all.min()) * 0.03
x_range = (x_all.min() - x_pad, x_all.max() + x_pad)
y_range = (y_all.min() - y_pad, y_all.max() + y_pad)

sc = None
for k, (ix, iy, x_lbl, y_lbl, title) in enumerate(resolved):
    ax = axes[k]
    sc = ax.scatter(all_dists[ix], all_dists[iy],
                    c=times, cmap="viridis",
                    s=args.marker_size, alpha=args.alpha,
                    linewidths=0, rasterized=True)
    ax.set_title(title, loc="left")
    ax.set_xlabel(f"{x_lbl} (Å)")
    ax.set_ylabel(f"{y_lbl} (Å)")
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    sns.despine(ax=ax)

# Shared time colorbar on the right
cbar = fig.colorbar(sc, ax=axes.tolist(), location="right",
                    shrink=0.85, pad=0.02)
cbar.set_label("Time (ns)")

out = f"{args.prefix}_distances_scatter.{args.format}"
fig.savefig(out, dpi=args.dpi)
print(f"Saved {out}")
plt.close(fig)