import pandas as pd
import plotly.express as px
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--dfl", nargs='+')
parser.add_argument("--dfd", nargs='+')
args = parser.parse_args()

datas_lipid = pd.concat([pd.read_csv(df) for df in args.dfl])
datas_dihe = pd.concat([pd.read_csv(df) for df in args.dfd])

data = pd.merge(datas_lipid, datas_dihe)
data['dihe'] = data['dihe'].abs()
print(data)

fig = px.box(data, x='interface', y='dihe', hover_name='system', facet_row='ligand_type', color='occupied')
fig.write_html('lipid_vs_gate_dihe.html')