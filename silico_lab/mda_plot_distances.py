#!/usr/bin/env python3
"""Plot atom-pair distances from a CSV produced by distances_calc.py.

Layout: rows = groups, columns = replicas. For each group, the time-series
and histogram of every replica sit in the same row so they can be compared
directly at a glance."""

import argparse
import csv
from collections import defaultdict
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ── CLI ──────────────────────────────────────────────────────────────
p = argparse.ArgumentParser(description="Plot atom-pair distances from CSV")
p.add_argument("-i", "--input", required=True, help="Input CSV (from distances_calc.py)")
p.add_argument("-w", "--window", type=int, default=50, help="Smoothing window in frames (default: 50, 0=off)")
p.add_argument("-b", "--bins", type=int, default=60, help="Histogram bin count (default: 60)")
p.add_argument("--mean-range", type=float, nargs=2, default=[100.0, 200.0],
               metavar=("START_NS", "END_NS"),
               help="Time window (ns) over which the mean distance is taken "
                    "for the summary table (default: 100 200)")
p.add_argument("-o", "--prefix", default="distances", help="Output file prefix")
args = p.parse_args()

# ── Read CSV ─────────────────────────────────────────────────────────
# Header: time_ns, replica, "group|label", "group|label", ...
with open(args.input) as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = list(reader)

columns = header[2:]
groups = [c.split("|", 1)[0] for c in columns]
labels = [c.split("|", 1)[1] for c in columns]

# Group rows by replica, preserving insertion order
rep_rows = defaultdict(list)
for row in rows:
    rep_rows[row[1]].append(row)

unique_replicas = list(rep_rows.keys())

# Per-replica numpy arrays
replica_data = {}
for rep, rrows in rep_rows.items():
    times = np.array([float(r[0]) for r in rrows])
    dists = np.array([[float(x) for x in r[2:]] for r in rrows]).T  # (n_pairs, n_frames)
    replica_data[rep] = (times, dists)

# Unique groups in order of first appearance → one subplot row each
unique_groups = []
for g in groups:
    if g not in unique_groups:
        unique_groups.append(g)
group_to_row = {g: i + 1 for i, g in enumerate(unique_groups)}

n_rows = len(unique_groups)
n_cols = len(unique_replicas)

print(f"Loaded {len(labels)} pairs, {n_rows} groups, {n_cols} replicas: {unique_replicas}")
for rep, (t, _) in replica_data.items():
    print(f"  {rep}: {len(t)} frames, {t[0]:.2f}–{t[-1]:.2f} ns")


# ── Time-series figure ──────────────────────────────────────────────
def smooth(y, window=50):
    """Running mean with proper edge handling.

    `np.convolve(..., mode='same')` implicitly pads with zeros, which pulls
    the smoothed line toward 0 over the first/last window/2 frames. Here we
    divide by the *actual* number of samples that fell inside the kernel at
    each position (computed by convolving a unit signal with the same
    kernel). At the edges the window effectively shrinks; in the interior
    it's identical to a plain running mean."""
    if window <= 1:
        return y
    kernel = np.ones(window)
    summed = np.convolve(y, kernel, mode="same")
    counts = np.convolve(np.ones_like(y, dtype=float), kernel, mode="same")
    return summed / counts


fig = make_subplots(
    rows=n_rows, cols=n_cols,
    shared_xaxes="columns",  # all groups within a replica share its time axis
    #shared_yaxes="rows",     # same group compared across replicas → same y scale
    shared_yaxes="all",
    vertical_spacing=0.06, horizontal_spacing=0.04,
    row_titles=unique_groups, column_titles=unique_replicas,
)

colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

for j, (label, group) in enumerate(zip(labels, groups)):
    row_idx = group_to_row[group]
    color = colors[j % len(colors)]
    for col_idx, rep in enumerate(unique_replicas, start=1):
        times, all_dists = replica_data[rep]
        ref_d = all_dists[j, 0]
        # TODO: toggle initial value substraction
        raw = all_dists[j] - ref_d
        #raw = all_dists[j]
        show_in_legend = (col_idx == 1)  # one legend entry per pair, in first replica column

        # Raw trace (thin, semi-transparent)
        fig.add_trace(
            go.Scattergl(
                x=times, y=raw, mode="lines", name=label,
                line=dict(width=0.5, color=color), opacity=0.3,
                showlegend=show_in_legend,
                legendgroup=label,
                hovertemplate=f"{label} ({rep})<br>t=%{{x:.2f}} ns<br>d=%{{y:.2f}} Å<extra></extra>",
            ),
            row=row_idx, col=col_idx,
        )
        # Smoothed trace (bold, on top)
        if args.window > 1:
            fig.add_trace(
                go.Scattergl(
                    x=times, y=smooth(raw, args.window), mode="lines",
                    name=f"{label} (avg {args.window}f)", showlegend=False,
                    line=dict(width=2, color=color),
                    legendgroup=label,
                    hovertemplate=f"{label} ({rep}) smoothed<br>t=%{{x:.2f}} ns<br>d=%{{y:.2f}} Å<extra></extra>",
                ),
                row=row_idx, col=col_idx,
            )

# Reference line at 4 Å — typical salt-bridge / strong polar contact threshold
for r in range(1, n_rows + 1):
    for c in range(1, n_cols + 1):
        fig.add_hline(y=1.5, line=dict(width=1, color="gray", dash="dot"),
                      row=r, col=c)
        fig.add_hline(y=-1.5, line=dict(width=1, color="gray", dash="dot"),
                      row=r, col=c)

# y-axis title only on leftmost column; x-axis title only on bottom row
for r in range(1, n_rows + 1):
    fig.update_yaxes(title_text="Distance (Å)", row=r, col=1)
for c in range(1, n_cols + 1):
    fig.update_xaxes(title_text="Time (ns)", row=n_rows, col=c)

fig.update_layout(
    height=max(900, 250 * n_rows),
    width=max(900, 450 * n_cols) + 280,  # extra width to host the right-side legend
    template="plotly_white",
    legend=dict(orientation="v", yanchor="top", y=1.0, xanchor="left", x=1.02,
                font=dict(size=9)),
    title="Atom-Pair Distances",
)

out_html = f"{args.prefix}_distances_plots.html"
fig.write_html(out_html, include_plotlyjs="cdn")
print(f"Saved {out_html}")
fig.show()

# ── Histogram figure ────────────────────────────────────────────────
# Distributions of the absolute distance for each pair, overlaid within each
# (group, replica) cell. Shared x within rows so the distance range is
# directly comparable across replicas of the same group.
fig_hist = make_subplots(
    rows=n_rows, cols=n_cols,
    shared_xaxes="rows", shared_yaxes="rows",
    vertical_spacing=0.08, horizontal_spacing=0.04,
    row_titles=unique_groups, column_titles=unique_replicas,
)

for j, (label, group) in enumerate(zip(labels, groups)):
    row_idx = group_to_row[group]
    color = colors[j % len(colors)]
    for col_idx, rep in enumerate(unique_replicas, start=1):
        _, all_dists = replica_data[rep]
        show_in_legend = (col_idx == 1)
        fig_hist.add_trace(
            go.Histogram(
                x=all_dists[j], name=label,
                marker=dict(color=color), opacity=0.5,
                nbinsx=args.bins,
                histnorm="probability density",
                showlegend=show_in_legend,
                legendgroup=label,
                hovertemplate=f"{label} ({rep})<br>d=%{{x:.2f}} Å<br>density=%{{y:.3f}}<extra></extra>",
            ),
            row=row_idx, col=col_idx,
        )

for r in range(1, n_rows + 1):
    fig_hist.update_yaxes(title_text="Density", row=r, col=1)
for c in range(1, n_cols + 1):
    fig_hist.update_xaxes(title_text="Distance (Å)", row=n_rows, col=c)

for r in range(1, n_rows + 1):
    for c in range(1, n_cols + 1):
        fig_hist.add_vline(x=4.0, line=dict(width=1, color="gray", dash="dot"),
                           row=r, col=c)

fig_hist.update_layout(
    barmode="overlay",
    height=max(900, 250 * n_rows),
    width=max(900, 450 * n_cols) + 280,
    template="plotly_white",
    legend=dict(orientation="v", yanchor="top", y=1.0, xanchor="left", x=1.02,
                font=dict(size=9)),
    title="Atom-Pair Distance Distributions",
)

out_html_hist = f"{args.prefix}_distances_histograms.html"
fig_hist.write_html(out_html_hist, include_plotlyjs="cdn")
print(f"Saved {out_html_hist}")
fig_hist.show()

# ── Summary table ────────────────────────────────────────────────────
# Per pair × replica:    Δ = ⟨d⟩(MEAN_START..MEAN_END ns)  −  d(t=0)
# Include pairs where |Δ| > THRESHOLD in at least one replica.
# Δ is signed: positive = the pair drifted apart on average; negative = it
# moved closer. BS_1 (chains A/B) and BS_2 (chains C/D) versions of the same
# interaction share a row (canonical = chain prefixes stripped).
THRESHOLD = 1.25    # Å — minimum |Δ| to include in table
MEAN_START, MEAN_END = args.mean_range   # ns; from CLI --mean-range

# Per-pair × per-replica signed shift; NaN if window has no frames for a replica
delta = np.full((len(labels), len(unique_replicas)), np.nan)
start_times = []
for ci, rep in enumerate(unique_replicas):
    times_arr, dists = replica_data[rep]
    start_times.append(times_arr[0])
    mask = (times_arr >= MEAN_START) & (times_arr <= MEAN_END)
    if not mask.any():
        print(f"  warning: {rep} has no frames in {MEAN_START}-{MEAN_END} ns window")
        continue
    delta[:, ci] = dists[:, mask].mean(axis=1) - dists[:, 0]

if any(abs(t) > 1e-3 for t in start_times):
    print(f"  note: 'd(t=0)' uses first frame in CSV; per-replica start times = "
          f"{', '.join(f'{t:.1f}' for t in start_times)} ns")


def site_of(label):
    """Classify by chain identifiers in the label: {A,B}→BS_1, {C,D}→BS_2,
    anything else → None (e.g. cross-site interface pairs)."""
    if " ↔ " not in label:
        return None
    chains = set()
    for side in label.split(" ↔ "):
        if len(side) >= 2 and side[1] == ":":
            chains.add(side[0])
    if chains and chains <= {"A", "B"}:
        return "BS_1"
    if chains and chains <= {"C", "D"}:
        return "BS_2"
    return None


def canonical(label):
    """Strip leading 'X:' chain prefix from each side of the pair label so
    BS_1 and BS_2 variants of the same interaction collapse to one key."""
    if " ↔ " not in label:
        return label
    return " ↔ ".join(
        s[2:] if len(s) >= 2 and s[1] == ":" else s
        for s in label.split(" ↔ ")
    )


# Group pairs by canonical → {site: original_index}; track insertion order
pair_by_canonical = {}
canonical_order = []
unclassified = []   # pairs whose chains don't fit BS_1 or BS_2 — not in table

for j, label in enumerate(labels):
    s = site_of(label)
    if s is None:
        unclassified.append(j)
        continue
    can = canonical(label)
    if can not in pair_by_canonical:
        pair_by_canonical[can] = {}
        canonical_order.append(can)
    pair_by_canonical[can][s] = j


def _canonical_abs_peak(can):
    """Largest |Δ| across all sites × replicas of a canonical pair."""
    m = 0.0
    for j in pair_by_canonical[can].values():
        x = np.abs(delta[j])
        if np.any(~np.isnan(x)):
            m = max(m, np.nanmax(x))
    return m


# Filter: keep canonicals where any (site, replica) had |Δ| > threshold
sig_canonicals = [
    can for can in canonical_order if _canonical_abs_peak(can) > THRESHOLD
]
print(f"\n{len(sig_canonicals)}/{len(pair_by_canonical)} canonical pair(s) "
      f"with |⟨d⟩({MEAN_START:.0f}–{MEAN_END:.0f} ns) − d(t=0)| > {THRESHOLD} Å in ≥1 replica")
if unclassified:
    print(f"  {len(unclassified)} pair(s) not classifiable as BS_1 / BS_2 "
          f"(chains outside {{A,B}} and {{C,D}}) — omitted from table:")
    for j in unclassified:
        print(f"    skip: {labels[j]}")

if len(sig_canonicals) > 0:
    # Sort by peak |Δ| across both sites, descending
    sig_canonicals.sort(key=lambda c: -_canonical_abs_peak(c))

    SITE_ORDER = ["BS_1", "BS_2"]
    DASH = "—"
    NAN  = "n/a"
    HIT  = "#ffe0e0"   # light red for cells over threshold (either sign)
    BG   = "white"

    def fmt(v):
        """Signed Δ with explicit +/- sign; (parens) if below threshold; n/a if NaN."""
        if np.isnan(v):
            return NAN
        s = f"{v:+.2f}"
        return s if abs(v) > THRESHOLD else f"({s})"

    # Build per-(site, replica) value and colour columns, in display order
    rep_columns = []
    rep_color_cols = []
    for site in SITE_ORDER:
        for ci in range(len(unique_replicas)):
            col_vals, col_cols = [], []
            for can in sig_canonicals:
                if site in pair_by_canonical[can]:
                    j = pair_by_canonical[can][site]
                    v = delta[j, ci]
                    col_vals.append(fmt(v))
                    col_cols.append(
                        HIT if (not np.isnan(v) and abs(v) > THRESHOLD) else BG
                    )
                else:
                    col_vals.append(DASH)
                    col_cols.append(BG)
            rep_columns.append(col_vals)
            rep_color_cols.append(col_cols)

    cell_values = [sig_canonicals] + rep_columns
    cell_colors = [[BG] * len(sig_canonicals)] + rep_color_cols
    header_values = ["Pair"] + [f"{s} {r}" for s in SITE_ORDER for r in unique_replicas]

    fig_table = go.Figure(data=[go.Table(
        columnwidth=[260] + [70] * (len(unique_replicas) * len(SITE_ORDER)),
        header=dict(
            values=header_values,
            fill_color="#d0d0d0",
            align="left",
            font=dict(size=12),
        ),
        cells=dict(
            values=cell_values,
            fill_color=cell_colors,
            align="left",
            font=dict(size=11, family="monospace"),
            height=22,
        ),
    )])

    fig_table.update_layout(
        title=(f"Pairs with |⟨d⟩({MEAN_START:.0f}–{MEAN_END:.0f} ns) − d(t=0)| > {THRESHOLD} Å "
               f"in ≥1 replica  ({len(sig_canonicals)} canonical pair(s)).  "
               f"BS_1 = chains A/B, BS_2 = chains C/D.  "
               f"Cells: signed Δ in Å; (parens) = below threshold; — = site absent; n/a = empty window."),
        height=max(300, 50 * len(sig_canonicals) + 120),
        width=max(1100, 260 + 80 * len(unique_replicas) * len(SITE_ORDER)),
        template="plotly_white",
    )

    out_html_table = f"{args.prefix}_distances_summary.html"
    fig_table.write_html(out_html_table, include_plotlyjs="cdn")
    print(f"Saved {out_html_table}")
    fig_table.show()