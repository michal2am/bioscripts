# podsumowuje zakres z wyników z localeval
# całość robi dla pierwszej sekwencji z .pir (druga to template)

import argparse
from mm_homology_sequencereader import Sequence
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns




class ScopeEvaluateLocal(Sequence):

    def __init__(self, pir_file, com_seq, sub_seq, sta_res, model_names, chains, energy_files, resi_span):

        self.models_names = model_names
        self.chains = chains

        Sequence.__init__(self, pir_file, com_seq, sub_seq, sta_res)

        self.set_value_csv(com_seq[0], energy_files[1], ['template'])
        self.set_value_csv(com_seq[0], energy_files[0], self.models_names)

        models = self.get_sequence(com_seq[0], num=True, skip=('residue', 'subunit', 'helix', 'strand', 'loop'))

        models_cut = models.loc[(models['position'] > 29) & (models['position'] < 61) & models['subunit'].isin(self.chains)]

        self.both_chains_sums = []
        self.both_chains_singles = []

        for ch in self.chains:

            sep_ch = models_cut.loc[(models['subunit'] == ch)]                                                          # select subunit and create a copy not to mess whole sequence
            singles = sep_ch.copy()

            singles.insert(4, 'mod_avg', singles.loc[:, self.models_names].mean(axis=1))
            sums = pd.DataFrame(singles.iloc[:, 3:].sum()).T
            sums.index = [ch]
            # sums.index = [ch+"_energy"]

            self.both_chains_sums.append(sums)
            self.both_chains_singles.append(singles)

            singles.to_csv(self.pir_file.strip('.pir') + '_cutene.csv', mode='a')
            sums.to_csv(self.pir_file.strip('.pir') + '_cutene.csv', mode='a')

    def plot_sum(self, relative=True):

        if relative:

            # drop average and drop template after relative
            both_chains_sums = pd.concat(self.both_chains_sums).drop(["mod_avg"], axis=1)
            both_chains_sums = both_chains_sums.iloc[:, :].div(both_chains_sums.template, axis=0)
            print("RELATIVE")
            print(both_chains_sums)

        else:

            # drop average and template
            both_chains_sums = pd.concat(self.both_chains_sums).drop(["template", "mod_avg"], axis=1)
            print("ABSOLUTE")
            print(both_chains_sums)

        # magic to go from two rows (chains) to row for each calculated DOPE
        # both_chains_sums["subunit"] = both_chains_sums.index                                                          # to avoid ugly subunit naming from scrip call arguments
        both_chains_sums["subunit"] = ['A1_1', 'A1_2']
        both_chains_sums["template"] = ['6HUP', '6HUP']
        both_chains_sums.index = [0, 1]
        both_chains_sums = both_chains_sums.melt(id_vars=['subunit', 'template'])

        print(both_chains_sums)

        if relative:

            # give sane column names and save realative DOPE sums
            both_chains_sums.rename(columns={'variable': 'model', 'value': 'relative DOPE'}, inplace=True)
            both_chains_sums.to_csv(self.pir_file.strip('.pir') + '_relative_sums.csv')

            # plot single box for both chains, relative DOPE sums

            sns.set_context("paper", font_scale=2.0)

            sns.set_style({'xtick.bottom': False,
                           'xtick.top': False,
                           'ytick.left': True,
                           'ytick.right': False})

            sns.swarmplot(y='relative DOPE', hue='subunit', x='template', dodge=True, data=both_chains_sums,
                          palette=sns.color_palette("muted"), size=7)

            ax = sns.boxplot(y='relative DOPE',  hue='subunit', x='template', dodge=True, data=both_chains_sums,
                             palette=['#BBBBBB', '#DDDDDD'], linewidth=2)

            handles, labels = ax.get_legend_handles_labels()
            l = plt.legend(handles[2:4], labels[2:4], loc=0, borderaxespad=0.)

            sns.despine(bottom=True, offset=10, trim=True)
            fig = plt.gcf()
            fig.set_size_inches(12, 7.5)
            fig.savefig('absolute_DOPE', dpi=300)

            plt.show()

            # plot barplot for each model (both chains grouped)
            ax2 = sns.barplot(x='model', y='relative DOPE', hue='subunit', data=both_chains_sums)
            ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90)

            # relative DOPE sums, but sorted
            both_chains_sums.sort_values(['subunit', 'relative DOPE', ], ascending=[True, False], inplace=True)
            both_chains_sums.to_csv(self.pir_file.strip('.pir') + '_relative_sums_sorted.csv')

            sns.despine()
            plt.show()

        else:

            # give sane column names and save  DOPE sums
            both_chains_sums.rename(columns={'variable': 'model', 'value': 'DOPE'}, inplace=True)
            both_chains_sums.to_csv(self.pir_file.strip('.pir') + '_absolute_sums.csv')

            # plot two boxes for each chain, absolute DOPE sums

            sns.set_context("paper", font_scale=2.0)

            sns.set_style({'xtick.bottom': False,
                           'xtick.top': False,
                           'ytick.left': True,
                           'ytick.right': False})

            sns.swarmplot(y='DOPE', hue='subunit', x='template', dodge=True, data=both_chains_sums,
                          palette=sns.color_palette("muted"), size=7)

            ax = sns.boxplot(y='DOPE',  hue='subunit', x='template', dodge=True, data=both_chains_sums,
                             palette=['#BBBBBB', '#DDDDDD'], linewidth=2)

            handles, labels = ax.get_legend_handles_labels()
            l = plt.legend(handles[2:4], labels[2:4], loc=0, borderaxespad=0.)

            sns.despine(bottom=True, offset=10, trim=True)
            fig = plt.gcf()
            fig.set_size_inches(12, 7.5)
            fig.savefig('absolute_DOPE', dpi=300)

            plt.show()



    def plot_singles(self):

        pass

        #print(self.both_chains_singles)
        #A = self.both_chains_singles[0]
        #B = self.both_chains_singles[1]

        #print(pd.merge(A, B, on=['residue', 'position']))
        #print(pd.merge(A, B, on=['residue', 'position']))
        #pd.merge(A, B, on=['residue', 'position']).to_csv('test.csv')
        #print(pd.concat([ A, B ], keys=['A','B'],axis=1))


        # barplot dla oddzielnych chainów

        #for chain in self.both_chains_singles:
        #    print(chain)
        #sorted_sums = sums.T.sort_values("energy").T
        #ax1 = sns.barplot(data=sorted_sums, palette=sns.color_palette("coolwarm", 102))
        #ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # ScopeEvaluateLocal

    parser.add_argument('-m', '--model_names', nargs='+')
    parser.add_argument('-c', '--chains', nargs='+')
    parser.add_argument('-e', '--energy_files', nargs='+')
    parser.add_argument('-s', '--resi_span', nargs='+')

    # Sequence
    parser.add_argument('-p',  '--pir_file')
    parser.add_argument('-cs', '--com_seq',   nargs='+')
    parser.add_argument('-ss', '--sub_seq',   nargs='+')
    parser.add_argument('-r',  '--sta_res',   nargs='+', type=int)

    args = parser.parse_args()

    evaluate = ScopeEvaluateLocal( args.pir_file, args.com_seq, args.sub_seq, args.sta_res, args.model_names,
                                   args.chains, args.energy_files, args.resi_span)

    analysis = 'full'

    if analysis == "full":

        evaluate.plot_sum(relative=False)
        evaluate.plot_sum(relative=True)

        # evaluate.plot_singles()