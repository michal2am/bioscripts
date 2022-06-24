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

# print(datasets)
# print(datasets[datasets['Occupancy'] > 15])


def res_plot_system(data):

    sns.set_style()
    sns.set_context("talk")
    g = sns.relplot(
        data=datasets[datasets['Occupancy'] > 0],
        x='Residue Num', y='Occupancy', hue='Residue Chain TypePos',
        row='Residue Chain Type',
        col='Lipid',
        # palette=sns.xkcd_palette(["pale red", 'windows blue']),
        # height=4, aspect=2,
        # common_norm=False,
        # facet_kws={'sharey':False},
        col_order=['POPC', 'DPSM', 'DPGS', 'PUPE', 'PAPS', 'PUPI', 'POP2', 'CHOL'],
        # hue_order=['gaba', 'bicuculline', 'flumazenil', 'phenobarbital', 'etomidate', 'propofol', 'diazepam', 'zolpidem']
    )

    g.set_titles("{row_name} | {col_name}")
    # g.set(xlim=(-1.5, 13), xlabel=r'distance to nearest lipid atom [$\AA$]', yticklabels=[], ylabel=None, yticks=[])
    # g.fig.suptitle(title)
    # g.fig.subplots_adjust(top=0.8)
    # g.despine(trim=False, left=True)
    # g.set_axis_labels("", "")
    # g.legend.set_title("")

    # plt.savefig('lipid_distances_sns_{}.png'.format(file), dpi=300)
    plt.savefig('lipid_occupancy.png', dpi=300)
    plt.show()

#res_plot_system(datasets)

def res_plot_lipid(data, lipid, systems, parameter):

    sns.set_style()
    sns.set_context("talk")
    data = data[(data['Lipid'] == lipid) & (data['System'].isin(systems))]
    print(data)

    for subunit in data['Residue Chain Type'].unique():
        data_subunit = data[data['Residue Chain Type'] == subunit]
        print(data_subunit.sort_values(parameter, ascending=False)[['Residue', 'Residue ID', 'Residue Chain Type','Residue Chain TypePos', 'System', parameter]].iloc[0:25, :])

    g = sns.relplot(
        data=data[data[parameter] > 0],
        x='Residue Num', y=parameter, hue='Residue Chain TypePos',
        row='Residue Chain Type',
        col='System',
        # palette=sns.xkcd_palette(["pale red", 'windows blue']),
        # height=4, aspect=2,
        # facet_kws={'sharey':False},
        col_order=systems,
        # hue_order=[]
    )

    g.set_titles("{row_name} | {col_name}")
    g.set(xlabel=r'Residue number')
    g.fig.suptitle(lipid)
    g.fig.subplots_adjust(top=0.9)
    g.despine(trim=False)
    # g.set_axis_labels("", "")
    g.legend.set_title('subunit_position')

    plt.savefig('{}_{}.png'.format(lipid, parameter.replace(" ", "")), dpi=300)
    # plt.show()

lipids = ['POPC', 'PUPE', 'CHOL', 'PAPS', 'PUPI', 'POP2', 'DPSM', 'DPGS']
#lipids = ['CHOL']
parameter = 'Residence Time'
for lipid in lipids:
    res_plot_lipid(datasets, lipid, ['6x3z', '6x3s'], parameter)