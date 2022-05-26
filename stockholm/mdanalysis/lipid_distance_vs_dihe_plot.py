import pandas as pd
import plotly.express as px
import argparse
import seaborn as sns
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
# python ~/repos/bioscripts/stockholm/mdanalysis/lipid_detector_ploter.py --df *lipids_distances.csv

parser.add_argument("--dfdist", nargs='+')
parser.add_argument("--dfdihe", nargs='+')


args = parser.parse_args()

dist_file = args.dfdist
dihe_file = args.dfdihe

lipid_data = pd.DataFrame()
dihe_data = pd.DataFrame()

for df in dist_file:

    single_data = pd.read_csv(df)
    lipid_data = pd.concat([lipid_data, single_data], ignore_index=True)

for df in dihe_file:

    single_data = pd.read_csv(df)
    dihe_data = pd.concat([dihe_data, single_data], ignore_index=True)

dihe_data['new_dihe'] = dihe_data['dihe'].apply(lambda x: x + 360 if x < 0 else x)
data = lipid_data.merge(dihe_data)
data = data[data.frame > 1600]

print(data)

def selective_plot_dihe(data, selected_lig, title, file):

    data = data[data.ligand_type.isin(selected_lig)]

    sns.set_style()
    sns.set_context("talk")
    g = sns.displot(
        data=data, kind="kde",
        x='new_dihe', y='distance', hue='ligand_type',
        row='ligand_state',
        # palette=sns.xkcd_palette(["pale red", 'windows blue']),
        height=4, aspect=2,
        # common_norm=False,
        # facet_kws={'sharey':False},
        row_order = ['holo', 'apo'],
        hue_order=['gaba', 'bicuculline', 'flumazenil', 'phenobarbital', 'etomidate', 'propofol', 'diazepam','zolpidem']
    )

    g.set(xlim=(-10, 370), xlabel=r'$\beta$2M286 dihedral angle [$\degree$]', ylabel=r'distance to nearest lipid atom [$\AA$]')
    g.set_titles("{row_name}")
    g.fig.suptitle(title)
    g.fig.subplots_adjust(top=0.8)
    g.despine(trim=False)
    g.legend.set_title("")

    plt.savefig('dist_vs_dihe_sns_{}.png'.format(file), dpi=300)
    plt.show()


selective_plot_dihe(data, ['gaba', 'bicuculline', 'flumazenil', 'phenobarbital'],
               r'$\beta$/$\alpha$ interface' +'\n' +'(not directly binding ligands)', 'BA_free')

selective_plot_dihe(data, ['etomidate', 'propofol', 'diazepam', 'zolpidem'],
               r'$\beta$/$\alpha$ interface' +'\n' +'(directly binding ligands)', 'BA_occu')