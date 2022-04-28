import pandas as pd
import plotly.express as px
import argparse

parser = argparse.ArgumentParser()
# python ~/repos/bioscripts/stockholm/mdanalysis/lipid_vs_dihe_plot.py --dfl lipids_5/diazepam_* -dfd diazepam_*

parser.add_argument("--dfl", nargs='+')
parser.add_argument("--dfd", nargs='+')
args = parser.parse_args()

datas_lipid = pd.concat([pd.read_csv(df) for df in args.dfl])
# datas_lipid['ligand_type'] = datas_lipid['ligand_type'].str[0:-3]
datas_dihe = pd.concat([pd.read_csv(df) for df in args.dfd])

print(datas_lipid)
print(datas_dihe)

#data = pd.merge(datas_lipid, datas_dihe)
data = datas_lipid.merge(datas_dihe)    # in this case merging works by default finding correct 'on' and 'how'

# cut first xxxx frames out
data = data[data.frame > 1600]

#data['dihe'] = data['dihe'].abs()
data['new_dihe'] = data['dihe'].apply(lambda x: x + 360 if x < 0 else x)
print(data)

data.ligand_type = pd.Categorical(data.ligand_type, categories=['bicuculline', 'gaba', 'etomidate', 'propofol', 'zolpidem', 'diazepam', 'phenobarbital'])
# data.ligand_type = pd.Categorical(data.ligand_type, categories=['bicuculline', 'etomidate', 'propofol', 'phenobarbital'])
data = data.sort_values('ligand_type')
#data['occupied'] = data['occupied'].astype(str)

# original for single box
# fig = px.box(data, x='interface', y='dihe', hover_name='system', facet_row='ligand_type', facet_col='ligand_state', color='occupied')
fig = px.histogram(data, x='new_dihe', facet_row='ligand_type', facet_col='ligand_state', color='occupied', color_discrete_map={0: '#6cd4c5', 1: '#e27c7c'})
fig.update_yaxes(matches=None)

for annotation in fig.layout.annotations:
    annotation.text = annotation.text.split("=")[1]

fig.write_html('lipid_vs_gate_dihe.html')
fig.write_image('lipid_vs_gate_dihe.png', width=800, height=1200, scale=2)


for ligand in data['ligand_type'].unique():
    if ligand == 'bicuculline':
        continue
    selected = ligand
    selected_data = data[(data['ligand_type'] == selected) & (data['ligand_state'] == 'apo')]
    fig2 = px.histogram(selected_data, x='new_dihe', facet_row='system', facet_col='interface', color='occupied',
                        title=ligand, color_discrete_map={0: '#6cd4c5', 1: '#e27c7c'})
    fig2.update_yaxes(matches=None)
    fig2.write_html('{}_detailed_lipid_vs_gate_dihe.html'.format(selected))
    fig2.write_image('{}_detailed_lipid_vs_gate_dihe.png'.format(selected), width=800, height=800, scale=2)