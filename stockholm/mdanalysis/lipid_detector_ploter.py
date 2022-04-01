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
fig = px.scatter(lipid_data, x='frame', y='occupied', color='interface', facet_col='system')
fig.write_html(data_file + '.html')