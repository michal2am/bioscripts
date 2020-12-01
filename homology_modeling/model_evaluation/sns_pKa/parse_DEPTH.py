import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
from functools import reduce
import plotly.graph_objects as go
import plotly.offline as pyo



class ReadPropkaLog:
    """
    prokpka reader
    """

    def __init__(self, file_names, pre_chains, post_chains):
        """
        propka log reader
        :param list file_names: file names
        :param list pre_chains: chains as in pdb/log file e.g. ['A', 'B' ... ]
        :param list post_chains: correct aliases for chains e.g ['A1_1', 'B2_2' ... ]
        """

        self.files = file_names
        self.chains_map = dict(zip(pre_chains, post_chains))
        self.results_human, self.results_long = self.read_files()

    def read_files(self):
        """
        reads propka messy log in row-wise manner
        :return: long and tidy format of pKa results
        """

        dfs = []

        for file in self.files:

            print("Reading {}".format(file))

            with open(file) as rf:

                parsed = []
                template, method, ligands = file[:-4].split('_')
                ligand = 'on' if ligands == 'ligand' else 'off'

                for row in rf:
                    if any(substring in row[0:4] for substring in ["ARG", "LYS", "ASP", "GLU", "HIS", "TYR", 'Ori', 'Swa']):
                        line = row.split()
                        if (len(line)) == 22:

                            if '*' in line[3]:
                                newline = [line[0], int(line[1]), self.chains_map[line[2]], np.NaN]
                            else:
                                newline = [line[0], int(line[1]), self.chains_map[line[2]], float(line[3])]

                            fullname = "{}{}{}".format(newline[0], newline[1], newline[2])
                            newline.insert(0, fullname)
                            newline.insert(-1, 'S')
                            newline.insert(3, newline[3][0:3])
                            newline[4] = newline[4][-3:]
                            newline.insert(-1, method)
                            newline.insert(-1, ligand)
                            parsed.append(newline)

                        elif (len(line)) == 23:

                            newline = [line[1][-3:], int(line[2].lstrip('O')), self.chains_map[line[3][0]], float(line[4][:-1])]
                            fullname = "{}{}{}".format(newline[0], newline[1], newline[2])
                            newline.insert(0, fullname)
                            newline.insert(-1, 'C1')
                            newline.insert(3, newline[3][0:3])
                            newline[4] = newline[4][-3:]
                            newline.insert(-1, method)
                            newline.insert(-1, ligand)
                            parsed.append(newline)

                            as_std = newline.copy()
                            as_std[5] = 'S'
                            parsed.append(as_std)

                        # TODO: line counting is dirty, change for Original/Swapped
                        #  and take into account triple sets present in Hibbs structures
                        #  and A1H/A1 naming scheme issue
                        ''''
                        elif (len(line)) == 24:

                            newline = [line[2][-3:], int(line[3].lstrip('O')), self.chains_map[line[4][0]], float(line[5][:-1])]    #! lstrip for PIO
                            fullname = "{}{}{}".format(newline[0], newline[1], newline[2])
                            newline.insert(0, fullname)
                            newline.insert(-1, 'C2')
                            newline.insert(3, newline[3][0:2])
                            newline[4] = newline[4][-3:]
                            newline.insert(-1, method)
                            newline.insert(-1, ligand)
                            parsed.append(newline)
                        # '''

                labels = ["fullname", "resname", "resnum", "chain", "chain_n", "type", "method", "ligand", template]
                df = pd.DataFrame.from_records(parsed, columns=labels)
                df.drop_duplicates(inplace=True)

            dfs.append(df)

        print("Reading files completed")

        results_human = reduce(lambda left, right: pd.merge(left, right, on=["fullname", "resname", "resnum", "chain", "chain_n", "type", "method", "ligand"], how='outer'), dfs)
        results_long = results_human.melt(id_vars=['fullname', 'resname', 'resnum', 'chain', 'chain_n', 'type', 'method', 'ligand']).copy()
        results_long.rename(columns={'variable': 'template', 'value': 'pKa'}, inplace=True)

        # print(results_human)
        # print(results_long)

        # results_long.to_csv('test.csv')

        print("Files parsed into data frame")

        return results_human, results_long


class ReadDepthLog:
    """
    depth reader
    """

    def __init__(self, file_names, pre_chains, post_chains):
        """
        depth log reader
        :param list file_names: file names
        :param list pre_chains: chains as in pdb/log file e.g. ['A', 'B' ... ]
        :param list post_chains: correct aliases for chains e.g ['A1_1', 'B2_2' ... ]
        """

        self.files = file_names
        self.chains_map = dict(zip(pre_chains, post_chains))
        self.results_human, self.results_long = self.read_files()

    def read_files(self):
        """
        reads depth log as whole csv
        :return: long and tidy format of pKa results
        """

        raws = []

        for file in self.files:

            template, method, ligands = file[:-4].split('_')
            ligand = 'on' if ligands == 'ligand' else 'off'

            raw_results = pd.read_csv(file, sep='\t')
            raw_results.drop(axis=1, labels=['MC_DEPTH', 'polar_SC_DEPTH', 'HB_N', 'HB_O', 'EE', 'total_SC_ASA_percent'], inplace=True)
            raw_results.insert(loc=0, column='method', value=method)
            raw_results.insert(loc=0, column='template', value=template)
            raw_results.insert(loc=0, column='type', value='S')
            raw_results.insert(loc=0, column='ligand', value=ligand)
            raw_results['fullname'] = raw_results.apply(lambda row: row.res_name + str(row.res_num) + self.chains_map[row.chain[-1]], axis=1)
            raw_results['chain'] = raw_results.apply(lambda row: row.fullname[-6:-4], axis=1)
            raw_results['chain_n'] = raw_results.apply(lambda row: row.fullname[-3:], axis=1)
            raw_results = raw_results.reindex(columns=['fullname', 'res_name', 'res_num', 'chain', 'chain_n', 'type', 'method', 'template', 'ligand', 'predicted_pKa'])
            raw_results = raw_results.rename(columns={'res_name': 'resname', 'res_num': 'resnum', 'predicted_pKa': 'pKa'})
            raws.append(raw_results)

        results_long = pd.concat(raws, axis=0)

        results_human = results_long.pivot_table(index=['fullname', 'resname', 'resnum', 'chain', 'chain_n', 'type', 'method', 'ligand'],columns='template', values='pKa')
        results_human.reset_index(drop=False, inplace=True)
        results_human.columns.name = None

        #print(results_long)
        #print(results_human)

        return results_human, results_long


class Analyze:
    """
    parsing and plotting
    """

    def __init__(self, longs):
        """
        analyzing all results method-independently
        :param list longs: long format dataframes
        """
        print("Concatenating all data")
        self.results_long = pd.concat(longs, axis=0, sort=False)
        # self.results_long.to_csv('all_results_long.csv')

    def parse(self, chains):
        """
        finds residues of pKa in range 5.5 - 8.5 in at least one structure
        :param str chain: chain name e.g. 'B3'
        :return: None
        """

        print("Selecting residues in 5.5 - 8.5 pKa range")
        selected = self.results_long[self.results_long['pKa'].between(5.5, 8.5)].copy()
        selected = selected.reindex(columns=['fullname', 'resnum', 'resname', 'chain', 'chain_n', 'type', 'method',
                                             'template', 'ligand', 'pKa'])
        selected = selected.sort_values(by='resnum')

        for chain in chains:

            print("Parsing chain {}".format(chain))

            resnums_selected = selected[selected['chain'] == chain]['resnum'].drop_duplicates(keep='first').sort_values().values
            print("\nResidues in range 5.5 to 8.5 (at least one measurement in range):")
            print(resnums_selected)
            print(selected[((selected['chain'] == chain) & (selected['resnum'].isin(resnums_selected)))])
            selected[((selected['chain'] == chain) & (selected['resnum'].isin(resnums_selected)))].to_csv('candidates_' + chain + '.csv')

            coupled_resnums_selected = selected[(selected['chain'] == chain) & (selected['type'] != 'S')]['resnum'].drop_duplicates(keep='first').sort_values().values
            print("\nCoupled residues in range 5.5 to 8.5 (at least one measurement in range and at least once coupled:")
            print(coupled_resnums_selected)
            print(selected[((selected['chain'] == chain) & (selected['resnum'].isin(coupled_resnums_selected)))])
            selected[((selected['chain'] == chain) & (selected['resnum'].isin(coupled_resnums_selected)))].to_csv('candidates_coupled' + chain + '.csv')

            coupled_resnums_all = self.results_long[(self.results_long['chain'] == chain) & (self.results_long['type'] != 'S')]['resnum'].drop_duplicates(keep='first').sort_values().values
            print("\nCoupled residues (at least once coupled:")
            print(coupled_resnums_all)

            def select_for_diff(chain, resnums):
                """
                prepares a df for two-column difference analysis
                :param str chain: chain physiological name, e.g. 'B3'
                :param list resnums: residue numbers
                :return: df of selected only residues with fullname removed
                """

                selected_residues = self.results_long[((self.results_long['chain'] == chain) & (self.results_long['resnum'].isin(resnums)) & (self.results_long['type'] == 'S'))].copy()
                selected_residues = selected_residues.reindex(columns=['fullname','resnum', 'resname', 'chain', 'chain_n', 'type', 'method', 'template', 'ligand', 'pKa'])
                selected_residues.drop(labels=['fullname'], axis=1, inplace=True)

                return selected_residues

            def pivot_for_diff(selected_residues, property, sort):
                """
                find residues with property dependant pKa differences
                :param df selected_residues: residues to check, prepared by select_for_diff
                :param str property: property to calculate diff, e.g. "chain_n"
                :param list sort: properties to sort by results
                :return:
                """

                new_columns = selected_residues.columns.values.tolist()
                new_columns.remove(property)
                new_columns.remove('pKa')

                selected_pivot = selected_residues.pivot_table(index=new_columns, columns=property, values='pKa')
                selected_pivot['diff'] = selected_pivot.apply(lambda row: abs(row[-2] - row[-1]), axis=1)

                selected_tresh = selected_pivot[selected_pivot['diff'] >= 0.50].sort_index(level=sort)

                print("\nDifferences >= 0.50 by {}".format(property))
                print(selected_tresh.reset_index(drop=False)['resnum'].drop_duplicates(keep='first').sort_values().values)
                print(selected_tresh)
                selected_tresh.to_csv('candidates_' + property + '_' + chain + '.csv')

            selected_residues = select_for_diff(chain, resnums_selected)
            pivot_for_diff(selected_residues, 'ligand', ['template', 'resnum'])
            pivot_for_diff(selected_residues, 'chain_n', ['resnum'])

    def plot_all(self, resids, chain, method):

        print("Plotting multi-plot")

        selected = self.results_long[(self.results_long["resnum"].isin(resids)) & (self.results_long["chain"].isin(chain)) & (self.results_long["type"] == 'S') & (self.results_long["method"] == method)]
        selected = selected.sort_values(by='resnum')

        # plotly
        #fig = px.scatter_matrix(selected, dimensions=["sepal_width", "sepal_length", "petal_width", "petal_length"],
        #                       color="species")


        '''

        plotly_plot = px.scatter(selected, x="ligand", y="pKa", color="template", facet_col="resnum",
                                 template='presentation', facet_col_wrap=20, title=chain)
        plotly_plot.write_html('plotly_plot_{}_ligand.html'.format(chain))
        plotly_plot.show()

        plotly_plot = px.scatter(selected, x="chain_n", y="pKa", color="template", facet_col="resnum",
                                 template='presentation', facet_col_wrap=20, title=chain)
        plotly_plot.write_html('plotly_plot_{}_chain_n.html'.format(chain))
        plotly_plot.show()

        plotly_plot = px.box(selected, x="chain_n", y="pKa", facet_col="resnum",
                                 template='presentation', facet_col_wrap=5, title=chain)
        plotly_plot.write_html('plotly_plot_box_{}_chain_n.html'.format(chain))
        plotly_plot.show()
        '''

        # seaborn

        sns.set_style()
        sns.set_context("paper")

        pivoted = pd.pivot_table(selected, index=['resnum', 'chain_n', 'template'], columns='ligand', values='pKa')
        pivoted.reset_index(inplace=True)
        print(pivoted)

        g = sns.catplot(kind='point', data=selected, col='chain_n', row='ligand', x='resnum', y='pKa',
                        hue='template', hue_order=['6i53', '6huo', '6hup', '6huk', '6huj', '6hug', '6x3z', '6x3x', '6x3u', '6x3t', '6x3v', '6x3w', '6x3s', '6x40'],
                        legend_out=True, height=2, aspect=1, palette=sns.mpl_palette("tab20b", 14))

        g._legend.remove()

        # g.set_titles("{col_name} {row_name}")
        # g.set_titles("{col_name}")
        g.axes[0, 0].axes.set_ylim(1.5, 8.5)
        g.axes[0, 0].axes.set_yticks(ticks=[2, 4, 6, 8])
        #g.axes[0, 0].yaxis.set_major_locator(plt.MaxNLocator(4))

        g.despine(trim=True)

        # g.set_axis_labels("", "")
        g.set_xticklabels(rotation=90)
        #plt.legend(bbox_to_anchor=(1.2, 0.1), loc=2, borderaxespad=0.)

        plt.tight_layout()
        plt.savefig('seaborn_plot_{}_{}_lignad.png'.format(chain[0], chain[1]), dpi=600)
        plt.show()
        '''
        '''

    def single_plot(self, resnums):

        print("Plotting single plot")

        for resnum in resnums:

            # plotly
            plotly_plot = px.scatter(self.results_long[self.results_long['resnum'] == resnum],
                                     x='ligand', y='pKa', facet_col='chain_n', facet_row='method', color='template',
                                     template='presentation', width=800, height=400, )
            plotly_plot.write_html('plotly_plot_' + str(resnum) + '.html')

            # seaborn
            sns.set_style()
            sns.set_context("talk")

            sns.set_context("paper")
            g = sns.PairGrid(self.results_long[self.results_long['resnum'] == resnum],
                             y_vars="pKa",
                             x_vars=['chain_n', 'ligand'],
                             hue='template',
                             height=2, aspect=1,
                             )

            g.map(sns.pointplot, scale=1., errwidth=4)
            # g.set(ylim=(0, 1))
            # sns.despine(fig=g.fig, left=True)
            g.add_legend()

            g.savefig('seaborn_plot_' + str(resnum) + '.png')




read_propka_ligand_Aris = ReadPropkaLog(['6i53_propka_ligand.log', '6huk_propka_ligand.log', '6hup_propka_ligand.log',
                                         '6huo_propka_ligand.log', '6huj_propka_ligand.log', '6hug_propka_ligand.log'],
                                        ["B", "A", "E", "D", "C", "G"],
                                        ["B3A_pri", "A1A_pri", "B3A_bis", "A1A_bis", "Y2A_pri", "G0_000"])


read_propka_free_Aris = ReadPropkaLog(['6i53_propka_free.log', '6huk_propka_free.log', '6hup_propka_free.log',
                                       '6huo_propka_free.log', '6huj_propka_free.log', '6hug_propka_free.log'],
                                      ["B", "A", "E", "D", "C", "G"],
                                      ["B3A_pri", "A1A_pri", "B3A_bis", "A1A_bis", "Y2A_pri", "G0_000"])


read_propka_ligand_Hibbs = ReadPropkaLog(['6x3s_propka_ligand.log', '6x3t_propka_ligand.log', '6x3u_propka_ligand.log',
                                    '6x3v_propka_ligand.log', '6x3w_propka_ligand.log', '6x3x_propka_ligand.log',
                                    '6x3z_propka_ligand.log', '6x40_propka_ligand.log'],
                                   ["A", "C", "E", "B", "D", "G", "M", "F", "H", "I", "L", "J", "K"],
                                   ["B2H_pri", "B2H_bis", "Y2H_pri", "A1H_pri", "A1H_bis", "G0_000", "M0_000", "F0_000", "H0_000", "I0_000", "L0_000", "J0_000", "K0_000",])

read_propka_free_Hibbs = ReadPropkaLog(['6x3s_propka_free.log', '6x3t_propka_free.log', '6x3u_propka_free.log',
                                  '6x3v_propka_free.log', '6x3w_propka_free.log', '6x3x_propka_free.log',
                                  '6x3z_propka_free.log', '6x40_propka_free.log'],
                                 ["A", "C", "E", "B", "D", "G", "M", "F", "H", "I", "L", "J", "K"],
                                 ["B2H_pri", "B2H_bis", "Y2H_pri", "A1H_pri", "A1H_bis", "G0_000", "M0_000", "F0_000", "H0_000", "I0_000", "L0_000", "J0_000", "K0_000",])


analyze = Analyze([read_propka_ligand_Aris.results_long, read_propka_free_Aris.results_long,
                   read_propka_ligand_Hibbs.results_long, read_propka_free_Hibbs.results_long,])
                   #read_depth_ligand.results_long, read_depth_free.results_long])

#analyze.single_plot([155, 267])

analyze.parse(['B3A', 'B2H', 'A1A', 'A1H'])

# analyze.plot_all([69, 84, 95, 101, 153, 155, 182, 191, 267, 270], 'B3A', 'propka')
# analyzfe.plot_all([112, 274], 'B3A', 'propka')

# analyze.plot_all([52, 84, 95, 119, 153, 155, 165, 182, 267, 270], 'B2H', 'propka')
# analyze.plot_all([112, 274, 279], 'B2H', 'propka')

# 191 out
analyze.plot_all([52, 69, 84, 95, 101, 119, 153, 155, 165, 182, 267, 270], ['B2H', 'B3A'], 'propka')


# analyze.plot_all([44, 56, 59, 63, 98, 110, 138, 142, 149, 151, 166, 170, 199, 216, 218], 'A1A', 'propka')
# analyze.plot_all([42, 105, 156, 303], 'A1A', 'propka')

# analyze.plot_all([44, 56, 59, 63, 110, 138, 142, 151, 166, 170, 199, 218, 234], 'A1H', 'propka')
# analyze.plot_all([98, 117, 156, 225, 279, 303], 'A1H', 'propka')

# 234 out
analyze.plot_all([44, 56, 59, 63, 98, 110, 138, 142, 149, 151, 166, 170, 199, 216, 218], ['A1H', 'A1A'], 'propka')


# 98, Aris low, Hibbs high

'''

read_propka_ligand_Hibbs = ReadPropkaLog(['6x3v_propka_ligand.log', ],
                                   ["A", "C", "E", "B", "D", "G", "M", "F", "H", "I", "L", "J", "K"],
                                   ["B2_pri", "B2_bis", "Y2_pri", "A1_pri", "A1_bis", "G0_000", "M0_000", "F0_000", "H0_000", "I0_000", "L0_000", "J0_000", "K0_000",])

read_propka_free_Hibbs = ReadPropkaLog(['6x3v_propka_free.log', ],
                                   ["A", "C", "E", "B", "D", "G", "M", "F", "H", "I", "L", "J", "K"],
                                   ["B2_pri", "B2_bis", "Y2_pri", "A1_pri", "A1_bis", "G0_000", "M0_000", "F0_000", "H0_000", "I0_000", "L0_000", "J0_000", "K0_000",])

analyze = Analyze([read_propka_ligand_Hibbs.results_long, read_propka_free_Hibbs.results_long])
                   #read_depth_ligand.results_long, read_depth_free.results_long])

analyze.parse(['A1'])
'''

### OLD ###


'''
read_depth_ligand = ReadDepthLog(['6i53_depth_ligand.log', '6huk_depth_ligand.log', '6hup_depth_ligand.log',
                                  '6huo_depth_ligand.log', '6huj_depth_ligand.log', '6hug_depth_ligand.log'],
                                 ["B", "A", "E", "D", "C", "G"],
                                 ["B3_pri", "A1_pri", "B3_bis", "A1_bis", "Y2_pri", "G0_000"])
read_depth_free = ReadDepthLog(['6i53_depth_free.log', '6huk_depth_free.log', '6hup_depth_free.log',
                                '6huo_depth_free.log', '6huj_depth_free.log', '6hug_depth_free.log'],
                               ["B", "A", "E", "D", "C", "G"],
                               ["B3_pri", "A1_pri", "B3_bis", "A1_bis", "Y2_pri", "G0_000"])
'''


'''

read_propka_ligand_Hibbs = ReadPropkaLog(['6x3s_propka_ligand.log', '6x3t_propka_ligand.log',],
                                   ["A", "C", "E", "B", "D", "G", "M", "F", "H", "I", "L", "J", "K"],
                                   ["B2_pri", "B2_bis", "Y2_pri", "A1_bis", "A1_bis", "G0_000", "M0_000", "F0_000", "H0_000", "I0_000", "L0_000", "J0_000", "K0_000",])

read_propka_free_Hibbs = ReadPropkaLog(['6x3s_propka_free.log', '6x3t_propka_free.log', ],
                                 ["A", "C", "E", "B", "D", "G", "M", "F", "H", "I", "L", "J", "K"],
                                 ["B2_pri", "B2_bis", "Y2_pri", "A1_bis", "A1_bis", "G0_000", "M0_000", "F0_000", "H0_000", "I0_000", "L0_000", "J0_000", "K0_000",])

analyze = Analyze([
                   read_propka_ligand_Hibbs.results_long, read_propka_free_Hibbs.results_long])
analyze.exp_plot()
analyze.parse('B2')

'''

'''
read_depth_ligand = ReadDepthLog(['6i53_depth_ligand.log', '6huk_depth_ligand.log', '6hup_depth_ligand.log', '6huo_depth_ligand.log', '6huj_depth_ligand.log', '6hug_depth_ligand.log'], ["B", "A", "E", "D", "C", "G"], ["B3_pri", "A1_pri", "B3_bis", "A1_bis", "Y2_pri", "G0_000"])
read_depth_free = ReadDepthLog(['6i53_depth_free.log', '6huk_depth_free.log', '6hup_depth_free.log', '6huo_depth_free.log', '6huj_depth_free.log', '6hug_depth_free.log'], ["B", "A", "E", "D", "C", "G"], ["B3_pri", "A1_pri", "B3_bis", "A1_bis", "Y2_pri", "G0_000"])

analyze = Analyze([read_depth_ligand.results_long, read_depth_free.results_long, read_propka_ligand.results_long, read_propka_free.results_long])
analyze.exp_plot()
analyze.parse('B3')
'''


#analyze = Analyze([read_propka_ligand.results_long, read_propka_free.results_long])



#analyze.plot_all([48, 69, 84, 95, 101, 112, 119, 146, 153, 155, 182, 191, 267, 270, 274, 282], 'B3', 'depth')
#analyze.plot_all([48, 69, 84, 95, 101, 112, 119, 146, 153, 155, 182, 191, 267, 270, 274, 282], 'B3', 'propka')


#analyze.parse('A1')
#analyze.plot_all([40, 42, 44, 52,  56, 59, 63, 98, 102, 105, 110, 138, 142, 149, 151, 156, 166, 170, 199, 216, 218, 287], 'A1', 'depth')
#analyze.plot_all([44, 56, 59, 63, 98, 102, 105, 110, 142, 149, 151, 156, 166, 199, 216, 218, 287], 'A1', 'propka')


#analyze.parse('Y2')
#analyze.plot_all([54, 56, 71, 75, 110, 148, 150, 156, 161, 163, 168, 178, 285, 289, 297], 'Y2', 'depth')
#analyze.plot_all([54, 56, 71, 75, 110, 148, 150, 156, 168, 178, 285, 289, 297], 'Y2', 'propka')

'''
#g = sns.pairplot(self.results_long[self.results_long['chain'] == 'B3'], hue='template')
print(self.results_long[self.results_long['resnum'] == 155])
g = sns.PairGrid(self.results_long[self.results_long['resnum'] == 155], y_vars='pKa', x_vars=['chain_n', 'method', 'ligand'], hue='template')
g.map(sns.pointplot)
g.map(sns.swarmplot)
'''

'''
g1 = sns.catplot(kind='point', data=self.results_long[self.results_long['resnum'] == 155], row='method',
                col='chain_n', x='ligand', y='pKa', hue='template', legend_out=True, margin_titles=True, size=3)
sns.despine(trim=True)

g2 = sns.catplot(kind='point', data=self.results_long[self.results_long['resnum'] == 155], row='method',
                col='ligand', x='chain_n', y='pKa', hue='template', legend_out=True, margin_titles=True, size=3)
sns.despine(trim=True)

#plt.show()

g1.savefig('test.png')
'''