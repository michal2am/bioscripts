#!/usr/bin/env python3
"""Publication-quality plots of atom-pair distances using matplotlib + seaborn.

Reads the same CSV format as distances_plot.py (group|label headers) and writes
two static figures: time-series and overlaid histograms. PDF by default;
--format png/svg also available.
"""

import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="Publication-quality distance plots from CSV")
p.add_argument("-i", "--input", required=True, help="Input CSV (from distances_calc.py)")
p.add_argument("-w", "--window", type=int, default=50, help="Smoothing window in frames (default: 50, 0=off)")
p.add_argument("-b", "--bins", type=int, default=60, help="Histogram bin count (default: 60)")
p.add_argument("-o", "--prefix", default="distances", help="Output file prefix")
p.add_argument("--format", default="pdf", choices=["pdf", "png", "svg"],
               help="Output format (default: pdf)")
p.add_argument("--dpi", type=int, default=300, help="DPI for raster output (default: 300)")
p.add_argument("--width", type=float, default=7.0, help="Figure width in inches (default: 7.0)")
p.add_argument("--panel-height", type=float, default=1.8, help="Per-panel height in inches (default: 1.8)")
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

# Unique groups in order of first appearance → one panel each
unique_groups = []
for g in groups:
    if g not in unique_groups:
        unique_groups.append(g)

# Pair indices grouped by panel
panel_pairs = {g: [] for g in unique_groups}
for j, group in enumerate(groups):
    panel_pairs[group].append(j)

print(f"Loaded {len(times)} frames, {len(labels)} pairs, {len(unique_groups)} groups: {unique_groups}")

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
    "legend.fontsize": 7,
    "pdf.fonttype": 42,   # keep text editable in PDF (TrueType, not Type-3)
    "ps.fonttype": 42,
})

# Distinct colours within each panel; colours restart per panel so each
# panel is self-contained visually.
palette = sns.color_palette("tab10")


def smooth(y, window=50):
    """Running average with edge handling."""
    kernel = np.ones(window) / window
    return np.convolve(y, kernel, mode="same")


# ── Time-series figure ──────────────────────────────────────────────
n_rows = len(unique_groups)
fig_ts, axes = plt.subplots(n_rows, 1,
                            figsize=(args.width, args.panel_height * n_rows),
                            sharex=True, constrained_layout=True)
if n_rows == 1:
    axes = [axes]

for row, group in enumerate(unique_groups):
    ax = axes[row]
    js = panel_pairs[group]
    for k, j in enumerate(js):
        label = labels[j]
        ref_d = all_dists[j, 0]
        # match the plotly script: trace is Δ from initial frame
        raw = all_dists[j] - ref_d
        color = palette[k % len(palette)]
        if args.window > 1:
            ax.plot(times, raw, color=color, lw=0.4, alpha=0.25)
            ax.plot(times, smooth(raw, args.window), color=color, lw=1.2, label=label)
        else:
            ax.plot(times, raw, color=color, lw=1.0, label=label)
    ax.set_title(group, loc="left")
    ax.set_ylabel("Distance (Å)")
    ax.legend(loc="best", frameon=False,
              ncol=2 if len(js) > 4 else 1)
    sns.despine(ax=ax)

axes[-1].set_xlabel("Time (ns)")
fig_ts.suptitle("Atom-Pair Distances", fontweight="bold")

out_ts = f"{args.prefix}_distances_timeseries.{args.format}"
fig_ts.savefig(out_ts, dpi=args.dpi)
print(f"Saved {out_ts}")
plt.close(fig_ts)

# ── Histogram figure ────────────────────────────────────────────────
fig_h, axes_h = plt.subplots(n_rows, 1,
                             figsize=(args.width, args.panel_height * n_rows),
                             constrained_layout=True)
if n_rows == 1:
    axes_h = [axes_h]

for row, group in enumerate(unique_groups):
    ax = axes_h[row]
    js = panel_pairs[group]
    for k, j in enumerate(js):
        label = labels[j]
        color = palette[k % len(palette)]
        sns.histplot(
            x=all_dists[j], bins=args.bins, stat="density",
            color=color, alpha=0.35, label=label,
            element="step", lw=1.2, ax=ax,
        )
    ax.set_title(group, loc="left")
    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("Density")
    ax.legend(loc="best", frameon=False,
              ncol=2 if len(js) > 4 else 1)
    sns.despine(ax=ax)

fig_h.suptitle("Atom-Pair Distance Distributions", fontweight="bold")

out_h = f"{args.prefix}_distances_histograms.{args.format}"
fig_h.savefig(out_h, dpi=args.dpi)
print(f"Saved {out_h}")
plt.close(fig_h)