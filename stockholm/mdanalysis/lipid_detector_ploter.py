import pandas as pd
import plotly.express as px
import argparse

parser = argparse.ArgumentParser()
# python ~/repos/bioscripts/stockholm/mdanalysis/lipid_detector_ploter.py --df *lipids.csv

parser.add_argument("--df", nargs='+') # 'etomidate_apo_lipids'

args = parser.parse_args()

data_file = args.df
lipid_data = pd.DataFrame()


for df in data_file:

    single_data = pd.read_csv(df)
    lipid_data = pd.concat([lipid_data, single_data])

print(lipid_data)
lipid_data = lipid_data[lipid_data.frame > 1600]


cumulative_data = lipid_data.groupby(['system', 'interface', 'ligand_type', 'ligand_state']).apply(lambda x: x['occupied'].sum()/len(x))
cumulative_data = cumulative_data.reset_index(name='occupied')

print(cumulative_data)

# cumulative_data.ligand_type = pd.Categorical(cumulative_data.ligand_type, categories=['etomidate', 'propofol', 'zolpidem4', 'diazepam4', 'phenobarbital', 'flumazenilnogaba', 'bicuculline', ],)
cumulative_data.ligand_type = pd.Categorical(cumulative_data.ligand_type, categories=['propofol5nr', 'etomidate5nr','diazepam5nr', 'zolpidem5nr', 'phenobarbital5nr', 'bicuculline5nr'])
# cumulative_data.ligand_type = pd.Categorical(cumulative_data.ligand_type, categories=['propofol5nr', 'propofol6nr'])

cumulative_data = cumulative_data.sort_values('ligand_type')

fig_sum = px.box(cumulative_data, x='ligand_type', y='occupied', points='all', hover_name='system',
                 hover_data=['ligand_type', 'ligand_state'], facet_col='interface',
                 category_orders={'interface': ['1st_beta/alpha', '2nd_beta/alpha', 'gamma/beta', 'alpha/gamma', 'gamma/beta']},
                 color='ligand_state')

fig_sum.write_html('lipid_contacts' + '_sum.html')

#fig_circle = px.scatter_polar(cumulative_data, r="occupied", theta="interface", color="ligand_type", symbol="ligand_state")
#fig_circle.write_html('lipid_contacts' + '_circle.html')