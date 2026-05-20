import argparse
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

atom_resi_names = {
    'dist_abuElectro': [
        "T202_GC_BS1", "T202_GC_BS2", "R67_GC_BS1", "R67_GC_BS2",
        "T130_GC_BS1", "T130_GC_BS2", "E155_GN_BS1", "E155_GN_BS2",
        "Y97_GN_BS1", "Y97_GN_BS2", "S156_GN_BS1", "S156_GN_BS2",
    ],
    'dist_b9b10': [
        "R207_K196_1_BS1", "R207_K196_1_BS2", "R207_K196_2_BS1", "R207_K196_2_BS2",
        "Y205_V198_1_BS1", "Y205_V198_1_BS2", "Y205_V198_2_BS1", "Y205_V198_2_BS2",
        "G203_F200_1_BS1", "G203_F200_1_BS2", "G203_F200_2_BS1", "G203_F200_2_BS2",
    ],
    'dist_capping': [
        "T202_Y205_BS1", "T202_Y205_BS2", "F200_F46_BS1", "F200_F46_BS2",
    ],
    'dist_r207': [
        "E153_R207_BS1", "E153_R207_BS2", "E155_R207_BS1", "E155_R207_BS2",
    ],
}


def get_columns(filename):
    for key, cols in atom_resi_names.items():
        if key in filename:
            return cols
    raise ValueError(f"No column definition found for '{filename}'. Known prefixes: {list(atom_resi_names.keys())}")


def plots(filename, distances):
    fig_line = px.line(distances)

    for trace in fig_line.data:
        y0 = distances[trace.name].iloc[0]
        fig_line.add_trace(go.Scatter(
            x=[distances.index[0], distances.index[-1]],
            y=[y0, y0],
            mode='lines',
            line=dict(color=trace.line.color, dash='dash'),
            name=f"{trace.name} t=0",
            legendgroup=trace.name,
            showlegend=False,
        ))

    colors = {t.name: t.line.color for t in fig_line.data if t.name in distances.columns}

    fig_line.write_html(filename + '.html')
    fig_line.show()

    fig_hist = go.Figure()
    for col in distances.columns:
        fig_hist.add_trace(go.Histogram(
            x=distances[col],
            name=col,
            marker_color=colors[col],
            opacity=0.7,
            xbins=dict(size=0.1),
        ))

    fig_hist.update_layout(barmode='overlay', xaxis_title='Distance', yaxis_title='Count')
    fig_hist.write_html(filename + '_hist.html')
    fig_hist.show()


parser = argparse.ArgumentParser(description='Plot distance data with line and histogram views.')
parser.add_argument('--filenames', nargs='+', help='Data file(s) without extension, e.g. dist_abuElectro_1us_sys4')
args = parser.parse_args()

for filename in args.filenames:
    cols = get_columns(filename)
    distances = pd.read_csv(filename + '.dat', header=None, sep=r'\s+', names=cols)
    plots(filename, distances)


