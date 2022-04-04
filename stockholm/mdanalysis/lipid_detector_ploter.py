import pandas as pd
import plotly.express as px
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--df") # 'etomidate_apo_lipids'
args = parser.parse_args()

data_file = args.df

lipid_data = pd.read_csv(data_file + '.csv')
#lipid_data = lipid_data.astype({'frame': float, 'occupied': float})

print(lipid_data)

fig_traj = px.scatter(lipid_data, x='frame', y='occupied', color='interface', facet_col='system')
fig_traj.write_html(data_file + '_traj.html')

num_frame = max(lipid_data['frame'])
cumulative_data = lipid_data.groupby(['system', 'interface'])['occupied'].sum()/num_frame
cumulative_data = cumulative_data.to_frame().reset_index()

print(cumulative_data)
#cumulative_data.reset_index(inplace=True)

fig_sum = px.scatter(cumulative_data, x='interface', y='occupied', color='system')
fig_sum.write_html(data_file + '_sum.html')
