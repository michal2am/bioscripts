import pandas as pd
import plotly.express as px
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--dfl", nargs='+')
parser.add_argument("--dfd", nargs='+')
args = parser.parse_args()

datas_lipid = pd.concat([pd.read_csv(df) for df in args.dfl])
datas_dihe = pd.concat([pd.read_csv(df) for df in args.dfd])

#data = pd.merge(datas_lipid, datas_dihe)
data = datas_lipid.merge(datas_dihe)    # in this case merging works by default finding correct 'on' and 'how'
#data['dihe'] = data['dihe'].abs()
data['new_dihe'] = data['dihe'].apply(lambda x: x + 360 if x < 0 else x)
print(data)

# original for single box
# fig = px.box(data, x='interface', y='dihe', hover_name='system', facet_row='ligand_type', facet_col='ligand_state', color='occupied')
fig = px.histogram(data, x='new_dihe', hover_name='system', facet_row='ligand_type', facet_col='ligand_state', color='occupied')
fig.update_yaxes(matches=None)

fig.write_html('lipid_vs_gate_dihe.html')