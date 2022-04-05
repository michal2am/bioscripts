import pandas as pd
import plotly.express as px
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--df", nargs='+') # 'etomidate_apo_lipids'
args = parser.parse_args()

data_file = args.df
lipid_data = pd.DataFrame()


for df in data_file:

    single_data = pd.read_csv(df + '.csv')
    lipid_data = pd.concat([lipid_data, single_data])

print(lipid_data)

fig_traj = px.scatter(lipid_data, x='frame', y='occupied', color='interface', facet_col='system')
#fig_traj.write_html(data_file + '_traj.html')
fig_traj.write_html('lipid_contacts' + '_traj.html')


# num_frame = max(lipid_data['frame'])
# cumulative_data = lipid_data.groupby(['system', 'interface', 'ligand_type', 'ligand_state'])['occupied'].sum()/num_frame

cumulative_data = lipid_data.groupby(['system', 'interface', 'ligand_type', 'ligand_state']).apply(lambda x: x['occupied'].sum()/len(x))

cumulative_data = cumulative_data.reset_index(name='occupied')

print(cumulative_data)

# fig_sum = px.box(cumulative_data, x='interface', y='occupied', hover_name='system', points='all', facet_col='ligand_state', facet_row='ligand_type')
fig_sum = px.box(cumulative_data, x='interface', y='occupied', hover_name='system', points='all', facet_col='ligand_state', color='ligand_type')

# fig_sum.write_html(data_file + '_sum.html')
fig_sum.write_html('lipid_contacts' + '_sum.html')
