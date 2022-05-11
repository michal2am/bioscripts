import pandas as pd
import plotly.express as px
import argparse
import seaborn as sns
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
# python ~/repos/bioscripts/stockholm/mdanalysis/lipid_detector_ploter.py --df *lipids_distances.csv

parser.add_argument("--df", nargs='+') # 'etomidate_apo_lipids'

args = parser.parse_args()

data_file = args.df
lipid_data = pd.DataFrame()


for df in data_file:

    single_data = pd.read_csv(df)
    lipid_data = pd.concat([lipid_data, single_data], ignore_index=True)

print(lipid_data)
lipid_data = lipid_data[lipid_data.frame > 1600]

#lipid_data.ligand_type = pd.Categorical(lipid_data.ligand_type, categories=['bicuculline', 'gaba', 'etomidate', 'propofol', 'zolpidem', 'diazepam', 'phenobarbital'])
#lipid_data = lipid_data.sort_values('ligand_type')


def selective_plot_distance(data, selected_lig, selected_int, title, file):

    data = data[data.ligand_type.isin(selected_lig) & data.interface.isin(selected_int)]

    sns.set_style()
    sns.set_context("talk")
    g = sns.displot(
        data=data, kind="kde",
        x='distance', hue='ligand_type',
        row='ligand_state',
        # palette=sns.xkcd_palette(["pale red", 'windows blue']),
        height=4, aspect=2,
        facet_kws={'sharey':False}
    )

    g.set_titles("{row_name}")
    g.fig.suptitle(title)
    g.fig.subplots_adjust(top=0.8)
    g.despine(trim=True)
    # g.set_axis_labels("", "")
    g.legend.set_title("")

    plt.savefig('lipid_distances_sns_{}.png'.format(file), dpi=300)
    plt.show()


selective_plot_distance(lipid_data, ['gaba', 'bicuculline', 'flumazenil', 'phenobarbital'],
                        ['1st_beta/alpha', '2nd_beta/alpha'],
                        r'$\beta$/$\alpha$ interface' +'\n' +'(not directly bonding ligands)', 'BA_free')

selective_plot_distance(lipid_data, ['etomidate', 'propofol', 'diazepam', 'zolpidem'],
                        ['1st_beta/alpha', '2nd_beta/alpha'],
                        r'$\beta$/$\alpha$ interface' +'\n' +'(directly bonding ligands)', 'BA_occu')

selective_plot_distance(lipid_data, ['gaba', 'bicuculline', 'flumazenil', 'etomidate', 'propofol', 'diazepam', 'zolpidem'],
                        ['alpha/beta'],
                        r'$\alpha$/$\beta$ interface' +'\n' +'(not directly bonding ligands)', 'AB_free')

selective_plot_distance(lipid_data, ['phenobarbital'], ['alpha/beta'],
                        r'$\alpha$/$\beta$ interface' +'\n' +'(directly bonding ligands)', 'AB_occu')

selective_plot_distance(lipid_data, ['gaba', 'bicuculline', 'flumazenil', 'phenobarbital', 'etomidate', 'propofol', 'diazepam', 'zolpidem'],
                        ['alpha/gamma'],
                        r'$\alpha$/$\gamma$ interface' +'\n' +'(not directly bonding ligands)', 'AG_free')

selective_plot_distance(lipid_data, ['gaba', 'bicuculline', 'flumazenil', 'etomidate', 'propofol', 'zolpidem'],
                        ['gamma/beta'], r'$\gamma$/$\beta$ interface' +'\n' +'(not directly bonding ligands)', 'GB_free')

selective_plot_distance(lipid_data, ['phenobarbital', 'diazepam'], ['gamma/beta'],
                        r'$\gamma$/$\beta$ interface' +'\n' +'(directly bonding ligands)', 'GB_occu')