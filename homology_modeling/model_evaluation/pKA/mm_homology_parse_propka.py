import pandas as pd
import numpy as np
from functools import reduce
import matplotlib.pyplot as plt
import seaborn as sns

class ReadPropkaLog():

    def __init__(self, file_names, pre_chains, post_chains):
        self.files = file_names
        self.chains_map = dict(zip(pre_chains, post_chains))
        self.results_human, self.results_long = self.read_files()

    def read_files(self):

        dfs = []

        for file in self.files:

            with open(file) as rf:
                parsed = []
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
                            newline.insert(3, newline[3][0:2])
                            newline[4] = newline[4][-1]
                            newline.insert(-1, "propka")
                            parsed.append(newline)

                        elif (len(line)) == 23:

                            newline = [line[1][-3:], int(line[2]), self.chains_map[line[3][0]], float(line[4][:-1])]
                            fullname = "{}{}{}".format(newline[0], newline[1], newline[2])
                            newline.insert(0, fullname)
                            newline.insert(-1, 'C1')
                            newline.insert(3, newline[3][0:2])
                            newline[4] = newline[4][-1]
                            newline.insert(-1, "propka")
                            parsed.append(newline)

                        elif (len(line)) == 24:

                            newline = [line[2][-3:], int(line[3]), self.chains_map[line[4][0]], float(line[5][:-1])]
                            fullname = "{}{}{}".format(newline[0], newline[1], newline[2])
                            newline.insert(0, fullname)
                            newline.insert(-1, 'C2')
                            newline.insert(3, newline[3][0:2])
                            newline[4] = newline[4][-1]
                            newline.insert(-1, "propka")
                            parsed.append(newline)

                labels = ["fullname", "resname", "resnum", "chain", "chain_n", "type", "method", file[2:6]]
                df = pd.DataFrame.from_records(parsed, columns=labels)

            dfs.append(df)

        # human readable, template-as-column
        results_human = reduce(lambda left, right: pd.merge(left, right, on=["fullname", "resname", "resnum", "chain", "chain_n", "type", "method"], how='outer'), dfs)
        print(results_human)
        # panda readable, single-value-row
        results_long = results_human.melt(id_vars=['fullname', 'resname', 'resnum', 'chain', 'chain_n', 'type', 'method']).copy()
        results_long.rename(columns={'variable': 'template', 'value': 'pKa'}, inplace=True)
        print(results_long)
        return results_human, results_long

    def parse(self, chain):

        selected = self.results_long[self.results_long['pKa'].between(5.5, 8.5)].copy()
        all_resnum_selected = selected[selected['chain'] == chain]['resnum'].drop_duplicates(keep='first')
        print("Standards in range 5.5 to 8.5:")
        print(all_resnum_selected.sort_values().values)
        coupled_resnum_selected = selected[(selected['chain'] == chain) & (selected['type'] != 'S')]['resnum'].drop_duplicates(keep='first')
        print("Couples in range 5.5 to 8.5:")
        print(coupled_resnum_selected.sort_values().values)

        selected = selected.pivot_table(index=['fullname', 'resname', 'resnum', 'chain', 'chain_n', 'type', 'method'],
                                        columns='template', values='pKa')

        selected.reset_index(drop=False, inplace=True)
        selected.columns.name = None

        print(selected[(selected['type'] == 'S') & (selected['chain'] == chain)].sort_values(['resnum', 'chain_n']))
        print(selected[(selected['type'] != 'S') & (selected['chain'] == chain)].sort_values(['resnum', 'chain_n']))


    def plot_standards(self, resids, chain):

        sns.set_style()
        sns.set_context("talk")

        no_resides = len(resids)

        #fig1, f1_axes = plt.subplots(ncols=int(no_resides/2), nrows=2,
        #                             figsize=(4 * no_resides/2, 8))
        fig1, f1_axes = plt.subplots(ncols=9, nrows=1,
                                     figsize=(18, 6), sharey=True)
        #flat_axes = [item for sublist in f1_axes for item in sublist]
        flat_axes = f1_axes

        for ax, res in zip(flat_axes, resids):

            single_res = self.results_long[(self.results_long["resnum"] == res) & (self.results_long["chain"] == chain)]

            #sns.violinplot(ax=ax, x="chain_n", y="pKa", data=single_res, palette='pastel', inner='quart')
            #sns.boxplot(ax=ax, x="chain_n", y="pKa", data=single_res, palette='pastel')
            sns.pointplot(ax=ax, x="chain_n", y="pKa", data=single_res, hue="template", dodge=True)



            old = ax.get_yticks()
            print(old)
            new = [old[0], old[-1]]
            print(new)
            ax.set_yticks([4, 5, 6, 7])

            ax.legend_.remove()
            ax.set(ylabel='', xlabel=single_res['fullname'].iloc[0][:-4])

            sns.despine(ax=ax, trim=True)

        plt.tight_layout()
        plt.show()

        fig1.savefig("test.png", dpi=300)

    def plot_couples(self, resids, chain):

        sns.set_style()
        sns.set_context("talk")

        fig1, f1_axes = plt.subplots(ncols=2, nrows=3,
                                     figsize=(6, 9))
        flat_axes = [item for sublist in f1_axes for item in sublist]
        types = ['S', 'S', 'C1', 'C1', 'C2', 'C2']

        for ax, res, tp in zip(flat_axes, resids*3, types):

            single_res = self.results_long[(self.results_long["resnum"] == res) & (self.results_long["chain"] == chain) & (self.results_long['type'] == tp)]
            sns.violinplot(ax=ax, x="chain_n", y="pKa", data=single_res, palette='pastel', inner='quart')
            sns.pointplot(ax=ax, x="chain_n", y="pKa", data=single_res, hue="template", dodge=True)

            ax.legend_.remove()
            ax.set(ylabel='', xlabel=single_res['fullname'].iloc[0][:-4] + ' ' + tp)

            sns.despine(ax=ax, trim=True)

        sns.despine(trim=True)
        plt.tight_layout()
        plt.show()


test = ReadPropkaLog(['d_6i53.log', 'd_6huk.log', 'd_6hup.log', 'd_6huo.log', 'd_6huj.log', 'd_6hug.log'], ["B", "A", "E", "D", "C"], ["B3_1", "A1_1", "B3_2", "A1_2", "Y2_1"])

hum_res = test.results_human
hum_res.to_csv('propka_pka_hum.csv')
long_res = test.results_long
long_res.to_csv('propka_pka_long.csv')

parsed = test.parse('B3')
parsed = test.parse('A1')

test.plot_standards([48, 69, 84, 95, 101, 146, 153, 155, 191], 'B3')
#test.plot_standards([56, 59, 63, 98, 138, 142, 149, 151, 166, 216, 218, 250], 'A1')
#test.plot_standards([54, 56, 71, 110, 148, 150, 156, 163, 168, 178, 297], 'Y2')


#test.plot_couples(parsed, [298, 424], 'B3')
#test.plot_couples([267, 270], 'B3')

