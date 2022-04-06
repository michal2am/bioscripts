import pandas as pd
import plotly.express as px
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--df", nargs='+') # 'etomidate_apo_lipids'
args = parser.parse_args()

data_file = args.df
data = pd.DataFrame()

for df in data_file:

    single_data = pd.read_csv(df)
    data = pd.concat([data, single_data])

data['dihe'] = data['dihe'].abs()
print(data)

fig = px.box(data, x='interface', y='dihe', hover_name='system', facet_col='ligand_state', facet_row='ligand_type', color='system')
fig.write_html('gate_resid_dihe' + '.html')