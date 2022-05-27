import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--datasets", nargs='+')
args = parser.parse_args()

datasets = pd.concat([pd.read_csv(dataset) for dataset in args.datasets])
datasets.reset_index(inplace=True, drop=True)
datasets.to_csv('all_dataset.csv')

print(datasets)
print(datasets[datasets['Occupancy'] > 15])


def res_plot(data):

    sns.set_style()
    sns.set_context("talk")
    g = sns.relplot(
        data=datasets[datasets['Occupancy'] > 0],
        x='Residue Num', y='Residence Time', hue='Lipid',
        row='Residue Chain',
        col='Lipid',
        # palette=sns.xkcd_palette(["pale red", 'windows blue']),
        # height=4, aspect=2,
        # common_norm=False,
        # facet_kws={'sharey':False},
        col_order=['POPC', 'POSM', 'DPG3', 'POPE', 'POPS', 'POP2', 'CHOL'],
        # hue_order=['gaba', 'bicuculline', 'flumazenil', 'phenobarbital', 'etomidate', 'propofol', 'diazepam', 'zolpidem']
    )

    # g.set_titles("{row_name}")
    # g.set(xlim=(-1.5, 13), xlabel=r'distance to nearest lipid atom [$\AA$]', yticklabels=[], ylabel=None, yticks=[])
    # g.fig.suptitle(title)
    # g.fig.subplots_adjust(top=0.8)
    # g.despine(trim=False, left=True)
    # g.set_axis_labels("", "")
    # g.legend.set_title("")

    # plt.savefig('lipid_distances_sns_{}.png'.format(file), dpi=300)
    plt.show()

res_plot(datasets)