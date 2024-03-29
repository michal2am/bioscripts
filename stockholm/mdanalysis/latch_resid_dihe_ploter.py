import pandas as pd
import plotly.express as px
import argparse
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
# python ~/repos/bioscripts/stockholm/mdanalysis/latch_resid_dihe_ploter.py --df *dihe.csv

parser.add_argument("--df", nargs='+')

args = parser.parse_args()

data_file = args.df
dihe_data = pd.DataFrame()


for df in data_file:

    single_data = pd.read_csv(df)
    dihe_data = pd.concat([dihe_data, single_data], ignore_index=True)

dihe_data = dihe_data[dihe_data.frame > 1600]

dihe_data['new_dihe'] = dihe_data['dihe'].apply(lambda x: x + 360 if x < 0 else x)
print(dihe_data)


def selective_plot_dihe(data, selected_lig, title, file):

    data = data[data.ligand_type.isin(selected_lig)]

    sns.set_style()
    sns.set_context("talk")
    g = sns.displot(
        data=data, kind="kde",
        x='new_dihe', hue='ligand_type',
        row='ligand_state',
        # palette=sns.xkcd_palette(["pale red", 'windows blue']),
        height=4, aspect=2,
        common_norm=False,
        # facet_kws={'sharey':False},
        row_order = ['holo', 'apo'],
        hue_order=['gaba', 'bicuculline', 'flumazenil', 'phenobarbital', 'etomidate', 'propofol', 'diazepam','zolpidem']
    )

    g.set(xlim=(-10, 370), xlabel=r'$\beta$2M286 dihedral angle [$\degree$]', yticklabels=[], ylabel=None, yticks=[])
    g.set_titles("{row_name}")
    g.fig.suptitle(title)
    g.fig.subplots_adjust(top=0.8)
    g.despine(trim=False, left=True)
    g.legend.set_title("")

    plt.savefig('latch_dihe_sns_{}.png'.format(file), dpi=300)
    plt.show()


selective_plot_dihe(dihe_data, ['gaba', 'bicuculline', 'flumazenil', 'phenobarbital'],
               r'$\beta$/$\alpha$ interface' +'\n' +'(not directly binding ligands)', 'BA_free')

selective_plot_dihe(dihe_data, ['etomidate', 'propofol', 'diazepam', 'zolpidem'],
               r'$\beta$/$\alpha$ interface' +'\n' +'(directly binding ligands)', 'BA_occu')

