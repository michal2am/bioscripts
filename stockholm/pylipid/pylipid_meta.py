import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import statistics
import ast

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--datasets", nargs='+')
args = parser.parse_args()

datasets = pd.concat([pd.read_csv(dataset) for dataset in args.datasets])
datasets.reset_index(inplace=True, drop=True)

# those are to fix dictionaries read from pickles
datasets['Durations Sys1'] = datasets['Durations Sys1'].apply(lambda values: ast.literal_eval(values))
datasets['Durations Sys2'] = datasets['Durations Sys2'].apply(lambda values: ast.literal_eval(values))
datasets['Durations Sys3'] = datasets['Durations Sys3'].apply(lambda values: ast.literal_eval(values))
datasets['Durations Sys4'] = datasets['Durations Sys4'].apply(lambda values: ast.literal_eval(values))
datasets['Durations All'] = datasets['Durations All'].apply(lambda values: ast.literal_eval(values))

datasets.to_csv('all_dataset.csv')
# print(datasets)


def res_plot_lipid(data, lipid, systems, parameter):

    sns.set_style()
    sns.set_context("talk")
    data = data[(data['Lipid'] == lipid) & (data['System'].isin(systems))]
    # print(data)

    for subunit in data['Residue Chain Type'].unique():
        data_subunit = data[data['Residue Chain Type'] == subunit]
        # print(data_subunit.sort_values(parameter, ascending=False)[['Residue', 'Residue ID', 'Residue Chain Type','Residue Chain TypePos', 'System', parameter]].iloc[0:25, :])

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

    plt.savefig('{}_{}_{}.png'.format(lipid, parameter.replace(" ", ""), ''.join([str(sys) for sys in systems])), dpi=300)
    # plt.show()


def event_times(system, lipid, time_tresh, n_tresh):
    """
    checks lipid interaction on the basis of raw event binding times and assigns results as beta parameter
    :param system: pdb code for pylipid dataset and pdb file name and directory
    :param lipid: CG residue name
    :param time_tresh: minimum interaction time
    :param n_tresh: minimum number of events of given time
    :return:
    """

    def scan_times(residue):
        # if ((max(residue['Durations Sys1']) > 2) & (max(residue['Durations Sys2']) > 2) & (max(residue['Durations Sys3']) > 2) & (max(residue['Durations Sys4']) > 2)):
        filter_times = [t for t in residue['Durations All'] if t > time_tresh]
        if len(filter_times) > n_tresh:
            residue_ids.append(residue['Residue ID'])

    residue_ids = []
    datasets[(datasets['Lipid'] == lipid) & (datasets['System'] == system)].apply(scan_times, axis=1)
    selected_byevents = datasets[
        (datasets['Lipid'] == lipid) & (datasets['System'] == system) & (datasets['Residue ID'].isin(residue_ids))]

    print(selected_byevents)

    pdb = open('{}_CG_b3/{}_CG_b3.pdb'.format(system, system), 'r')
    pdb_lines = pdb.readlines()
    new_pdb_lines = []

    for line in pdb_lines:
        if 'PRO' in line:
            atom = line.strip().split()
            resname = atom[3]
            resid = int(atom[5])
            segname = atom[11]
            detected = selected_byevents[
                (selected_byevents['Residue Name'] == resname) & (selected_byevents['Residue Num'] == resid) & (
                            selected_byevents['Residue Chain'] == segname)]
            if len(detected) > 0:
                new_line = line.strip()[0:62] + '1.00' + line.strip()[66:] + '\n'
                new_pdb_lines.append(new_line)
            else:
                new_pdb_lines.append(line)
        else:
            new_pdb_lines.append(line)

    new_pdb = open('{}_event_times.pdb'.format(system), 'w')
    new_pdb.writelines(new_pdb_lines)
    new_pdb.close()


make_plots = False
print_tops = False
map_times = True

if make_plots:
    #lipids = ['POPC', 'PUPE', 'CHOL', 'PAPS', 'PUPI', 'POP2', 'DPSM', 'DPGS']
    lipids = ['CHOL']
    parameters = ['Occupancy', 'Residence Time']
    for parameter in parameters:
        for lipid in lipids:
            res_plot_lipid(datasets, lipid, ['6x3z', '6x3s'], parameter)
            res_plot_lipid(datasets, lipid, ['7qn7', '7qn9'], parameter)

if print_tops:
    selected = datasets[(datasets['Lipid'] == 'CHOL') & (datasets['System'] == '6x3z')]
    #selected = selected[selected['Binding Site ID'].isin([9])]
    selected = selected.sort_values(by='Residence Time', ascending=False).iloc[0:50, :]
    print(selected[['Residue', 'Residue ID', 'Residue Chain', 'Binding Site ID', 'Residence Time']])

if map_times:
    event_times('6x3z', 'CHOL', 2, 4)
