import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import statistics
import ast
import time


def res_plot_lipid(datasets, lipid, systems, parameter):
    """

    :param datasets:
    :param lipid:
    :param systems:
    :param parameter:
    :return:
    """

    selected_dataset = datasets[(datasets['Lipid'] == lipid) & (datasets['System'].isin(systems))]

    for system in systems:
        data_system = selected_dataset[selected_dataset['System'] == system]
        print('\nTop residues by {} for {} lipid in {} system:\n'.format(parameter, lipid, system))
        print(data_system.sort_values(parameter, ascending=False)[['Residue', 'Residue ID', 'Residue Chain Type','Residue Chain TypePos', 'System', 'Event Threshold', parameter]].iloc[0:25, :])

    g = sns.relplot(
        data=selected_dataset[selected_dataset[parameter] > 0],
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
    g.legend.set_title('Subunit type/position')

    plt.savefig('{}_{}_{}.png'.format(lipid, parameter.replace(" ", ""), ''.join([str(sys) for sys in systems])), dpi=300)
    # plt.show()


def event_times(datasets, system, lipid, time_thresh, n_thresh):
    """
    checks lipid interaction on the basis of raw event binding times and assigns results as beta parameter and boolean threshold in dataset
    :param datasets: universal df with data
    :param system: pdb code for pylipid dataset and pdb file name and directory
    :param lipid: CG residue name
    :param time_thresh: minimum interaction time
    :param n_thresh: minimum number of events of given time
    :return:
    """

    def scan_times(residue):
        # if ((max(residue['Durations Sys1']) > 2) & (max(residue['Durations Sys2']) > 2) & (max(residue['Durations Sys3']) > 2) & (max(residue['Durations Sys4']) > 2)):
        filter_times = [t for t in residue['Durations All'] if t >= time_thresh]
        if len(filter_times) >= n_thresh:
            if residue['Residence Time'] >= 3:
                # print(residue)
                residue_ids.append(residue['Residue ID'])

    def map_pdb():
        pdb = open('../{}_CG_b3/{}_CG_b3.pdb'.format(system, system), 'r')
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

        new_pdb = open('{}_{}_event_times.pdb'.format(lipid, system), 'w')
        new_pdb.writelines(new_pdb_lines)
        new_pdb.close()

    residue_ids = []
    datasets[(datasets['Lipid'] == lipid) & (datasets['System'] == system)].apply(scan_times, axis=1)
    datasets.loc[((datasets['Lipid'] == lipid) & (datasets['System'] == system)),'Event Threshold'] = 0
    datasets.loc[((datasets['Lipid'] == lipid) & (datasets['System'] == system) & (datasets['Residue ID'].isin(residue_ids))),'Event Threshold'] = 1

    selected_byevents = datasets[((datasets['Event Threshold'] == 1) & (datasets['Lipid'] == lipid) & (datasets['System'] == system))]

    print('\n Residues above threshold for {} lipid in {} system:\n'.format(lipid, system))
    print(selected_byevents[['Residue', 'Residue ID', 'Residue Chain TypePos', 'System', 'Event Threshold']])

    map_pdb()


def cumulative_threshold_plot(datasets, lipids):
    """
    prints results of threshold analysis
    :param datasets: universal df with data
    :param lipids: plots and prints done only for selected
    :return:
    """

    # by system

    cumulative_system = datasets.groupby(['System', 'State', 'Lipid'], as_index=False)['Event Threshold'].sum()

    print('\nResidues above threshold, for each subunit system/state/lipid:\n')
    print(cumulative_system[cumulative_system['Lipid'].isin(lipids)])

    g = sns.catplot(
        data=cumulative_system[cumulative_system['Lipid'].isin(lipids)],
        x='Lipid', y='Event Threshold', hue='System',
        kind='bar',
        col='State',
        height=4, aspect=1,
        hue_order=['6x3z', '6x3s', '7qn7', '7qn9', '7qn8', '7qna', '7qnb', '6a96']
    )

    g.set_titles("{col_name}")
    g.despine(trim=False)
    g.set_axis_labels('Lipid', '# of contacting residues')
    plt.savefig('lipid_threshold_counts_system.png', dpi=300)

    # by subunit

    # manual mean calculation
    cumulative_subunit = datasets.groupby(['Residue Chain TypePos','Residue Chain Type', 'System', 'State', 'Lipid'], as_index=False)['Event Threshold'].sum()
    cumulative_subunit = cumulative_subunit.groupby(['Residue Chain Type', 'State', 'Lipid'], as_index=False)['Event Threshold'].mean()

    # manual n of subunits calculation
    subunits_no = datasets.groupby(['Residue Chain TypePos', 'Residue Chain Type', 'State', 'System'], as_index=False)['Event Threshold'].count()
    subunits_no = subunits_no.groupby(['Residue Chain Type', 'State'], as_index=False)['Event Threshold'].count()
    subunits_no.rename(columns={'Event Threshold': 'N subunits'}, inplace=True)

    cumulative_subunit_std = cumulative_subunit.merge(subunits_no)
    print('\nResidues above threshold, for each subunit type/state/lipid:\n')
    print(cumulative_subunit_std[cumulative_subunit_std['Lipid'].isin(lipids)])

    # prepare for means using catplot
    cumulative_subunit_tomean = datasets.groupby(['Residue Chain TypePos', 'Residue Chain Type', 'System', 'State', 'Lipid'], as_index=False)['Event Threshold'].sum()

    g = sns.catplot(
        data=cumulative_subunit_tomean[cumulative_subunit_tomean['Lipid'].isin(lipids)],
        x='Lipid', y='Event Threshold', hue='Residue Chain Type',
        kind='bar',
        col='State',
        ci=68,
        height=4, aspect=1,
        # hue_order=['POPC', 'PUPE', 'CHOL', 'PAPS', 'PUPI', 'POP2', 'DPSM', 'DPGS']
    )

    g.set_titles("{col_name}")
    g.legend.set_title('Subunit type')
    g.despine(trim=False)
    g.set_axis_labels('Lipid', '# of contacting residues')
    plt.savefig('lipid_threshold_counts_subunit.png', dpi=300)

    # by subunit position

    cumulative_subunit_position = datasets.groupby(['Residue Chain TypePos', 'System', 'State', 'Lipid'], as_index=False)['Event Threshold'].sum()
    print('\nResidues above threshold, for each subunit/position/system/state/lipid:\n')
    print(cumulative_subunit_position[cumulative_subunit_position['Lipid'].isin(lipids)])

    g = sns.catplot(
        data=cumulative_subunit_position[cumulative_subunit_position['Lipid'].isin(lipids)],
        y='Event Threshold', x='Residue Chain TypePos', hue='System',
        kind='bar',
        row='State',
        col='Lipid',
        height=4, aspect=1,
        hue_order=['6x3z', '6x3s', '7qn7', '7qn9', '7qn8', '7qna', '7qnb', '6a96']
    )

    g.set_titles("{row_name} | {col_name}")
    g.set_xticklabels(rotation=90)
    g.despine(trim=False)
    g.set_axis_labels('Subunit type/position', '# of contacting residues')
    plt.savefig('lipid_threshold_counts_subunit_position.png', dpi=300)


parser = argparse.ArgumentParser()
parser.add_argument('--parse_datasets')#, action=argparse.BooleanOptionalAction)
parser.add_argument('-d', '--datasets', nargs='+')
parser.add_argument('--map_thresholds')#, action=argparse.BooleanOptionalAction)

parser.add_argument('-systems', nargs='+', default=['6x3z', '6x3s', '7qn7', '7qn9', '7qn8', '7qna', '7qnb', '6a96'])
parser.add_argument('-lipids', nargs='+', default=['POPC', 'PUPE', 'CHOL', 'PAPS', 'PUPI', 'POP2', 'DPSM', 'DPGS'])
args = parser.parse_args()


if args.parse_datasets:
    # works on each csv, not dependent on systems/lipids parameters
    # read all datafiles from pylipid
    datasets = pd.concat([pd.read_csv(dataset) for dataset in args.datasets])
    datasets.reset_index(inplace=True, drop=True)

    # adds some general data
    states = {'6x3z': 'desensitized', '6x3s': 'closed', '7qn7': 'desensitized', '7qn9': 'closed',
              '7qn8': 'desensitized', '7qna': 'closed', '7qnb': 'closed', '6a96': 'open'}
    datasets['State'] = datasets['System'].apply(lambda sys: states[sys])
    datasets['Chain TypePosSys'] = datasets['Residue Chain TypePos'] + datasets['System']

    assemblies = {'6x3z': 'b2a1b2a1g2', '6x3s': 'b2a1b2a1g2', '7qn7': 'a4b3b3b3d', '7qn9': 'a4b3b3b3d',
                  '7qn8': 'b3b3b3b3d', '7qna': 'a4b3g2b3b3', '7qnb': 'g2b3g2b3b3', '6a96': 'a5b3b3b3b3'}
    datasets['Assembly'] = datasets['System'].apply(lambda sys: assemblies[sys])

    datasets.to_csv('all_dataset.csv', index=False)
else:
    datasets = pd.read_csv('all_dataset.csv')

if args.map_thresholds:
    # works only on selected systems/lipids
    # those are to fix dictionaries read from pickles
    datasets['Durations Sys1'] = datasets['Durations Sys1'].apply(lambda values: ast.literal_eval(values))
    datasets['Durations Sys2'] = datasets['Durations Sys2'].apply(lambda values: ast.literal_eval(values))
    datasets['Durations Sys3'] = datasets['Durations Sys3'].apply(lambda values: ast.literal_eval(values))
    datasets['Durations Sys4'] = datasets['Durations Sys4'].apply(lambda values: ast.literal_eval(values))
    datasets['Durations All'] = datasets['Durations All'].apply(lambda values: ast.literal_eval(values))

    # adds Event Threshold metric and maps pdbs betas
    for lipid in args.lipids:
        for system in args.systems:
            event_times(datasets, system, lipid, time_thresh=3, n_thresh=2)

    datasets.to_csv('all_dataset.csv', index=False)

sns.set_style()
sns.set_context("talk")
sns.set_palette('muted')


print_tops = True

if print_tops:
    selected = datasets[(datasets['Lipid'] == 'CHOL') & (datasets['System'] == '6x3z') & (datasets['Residue Name'].isin(['ARG', 'LYS']))]
    # selected = selected[selected['Binding Site ID'].isin([9])]
    # selected = selected.sort_values(by='Residence Time', ascending=False).iloc[0:50, :]
    selected = selected.sort_values(by='Occupancy', ascending=False).iloc[0:50, :]
    print(selected[['Residue', 'Residue ID', 'Residue Chain Type', 'Residue Chain TypePos', 'System', 'Occupancy', 'Residence Time', 'Event Threshold']])

# works on selected lipids and systems
parameters = ['Occupancy', 'Residence Time']
for parameter in parameters:
    for lipid in args.lipids:
        # some magic to group system by assemblies in correct order
        for assembly in datasets.loc[datasets['System'].isin(args.systems), 'Assembly'].unique():
            same_assembly_systems = datasets.loc[(datasets['System'].isin(args.systems)) &
                                                 (datasets['Assembly'] == assembly), 'System'].unique()
            ordered_systems = [sys for sys in args.systems if sys in same_assembly_systems]
            res_plot_lipid(datasets, lipid, ordered_systems, parameter)

# works on selected lipids, all systems
cumulative_threshold_plot(datasets, args.lipids)

