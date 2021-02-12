# meta analysis of WT in SPBUN

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import researchpy as rp


class MetaWT:

    def __init__(self):

        self.data = pd.read_csv("meta_wt_raw.csv")
        self.data = self.data.loc[:, ~self.data.columns.str.contains('^Unnamed')]
        self.data = self.data.drop(columns=['ds_s', 'da_s', 'dal_s'])

        parsed = self.data.pivot_table(index=['who', 'file '], columns=["type"]).copy()
        self.column_order = ['RT', 'RTln', 'FR10', 'FR300', 'FR500', 'ds_t2', 'ds_A2', 'ds_t1', 'ds_A1', 'ds_c', 'ds1_t', 'ds1_c', 'dal_t2', 'dal_A2', 'dal_t1', 'dal_A1', 'dal_m', 'da_t2', 'da_A2', 'da_t1', 'da_A1', 'da_m']
        parsed = parsed.reindex_axis(self.column_order, axis=1, level=0)

        print(parsed)
        print(parsed.describe())
        parsed.describe().to_csv('summary.csv')

        print(parsed.loc["Marek", :])
        print(parsed.loc["Marek", :].describe())
        parsed.loc["Marek", :].describe().to_csv('marek_summary.csv')

    def exp_lc_test(self):

        lcs = self.data[(self.data['type'] == 'lc')]
        lcs.reset_index(inplace=True)
        exps = self.data[(self.data['type'] == 'exp')]
        exps.reset_index(inplace=True)

        p_values = {}
        for param in lcs.columns[4:]:
            #print(param)
            descriptives, results = rp.ttest(lcs[param].dropna(how='all'), exps[param].dropna(how='all'))
            #print(descriptives)
            #print(results)
            p_values[param] = [results.iloc[3, 1], descriptives.loc[0, 'Mean'], descriptives.loc[1, 'Mean']]

        statistics = pd.DataFrame(p_values, index=['p-value', 'LC', 'EXP'])
        statistics = statistics.reindex_axis(self.column_order, axis=1)
        print(statistics)
        statistics.to_csv('statistics.csv')

        # print(stats.levene(lcs['FR10'].dropna(how='all'), exps['FR10'].dropna(how='all')))
        # print(stats.ttest_ind(lcs['FR10'].dropna(how='all'), exps['FR10'].dropna(how='all')))

    def plot_violin(self, param, yname):

        fig = plt.figure(figsize=(12, 5))
        grid = plt.GridSpec(1, 5, wspace=0.5, hspace=0.)
        ax1 = fig.add_subplot(grid[0, 0:4])
        ax2 = fig.add_subplot(grid[0, 4], sharey=ax1)

        g1 = sns.swarmplot(ax=ax1, x="who", y=param, data=self.data, hue="type", palette='deep', size=8)
        g1 = sns.violinplot(ax=ax1, x="who", y=param, data=self.data, hue="type", inner="quart", split=True,
                            palette='pastel')
        g1.legend_.remove()
        g1.set(ylabel=yname, xlabel='')

        g2 = sns.swarmplot(ax=ax2, y=self.data[param], x=[""] * len(self.data[param]), hue=self.data["type"],           # dirty hack without data='smth' to have hue with single 'x'
                           dodge=False, palette='deep', size=8)
        g2 = sns.violinplot(ax=ax2, y=self.data[param], x=[""] * len(self.data[param]), hue=self.data["type"],
                            inner='quart', split=True, palette='pastel', size=8)
        g2.legend_.remove()
        g2.set(ylabel='')
        g2.set(xticks=[])

        sns.despine(ax=ax1, offset=5, trim=True)
        sns.despine(ax=ax2, offset=5, trim=True, bottom=True, left=False)

        #plt.tight_layout()
        plt.show()

        fig.savefig(yname+"_all.png", dpi=300)

    def plot_ds(self):

        for param, labe in zip(['ds_A1', 'ds_A2', 'ds_c', 'ds_t1', 'ds_t2'], ['desensitization A slow',
                                                                              'desenitization A fast',
                                                                              'desensitization c',
                                                                              'desensitization t slow',
                                                                              'desensitization t fast']):
            self.plot_violin(param, labe)

    def plot_ds1(self):

        for param, labe in zip(['ds1_t', 'ds1_c'], ['desensitization (1 component) t', 'desensitization (1 component) c']):
            self.plot_violin(param, labe)

    def plot_fr(self):

        for param in ['FR10', 'FR300', 'FR500']:
            self.plot_violin(param, param)

    def plot_rt(self):

        for param in ['RT', 'RTln']:
            self.plot_violin(param, param)

    def plot_dea_short(self):

        for param, labe in zip(['da_A1', 'da_A2', 'da_t1', 'da_t2', 'da_m'], ['deactivation (short pulse) A slow',
                                                                              'deactivation (short pulse) A fast',
                                                                              'deactivation (short pulse) t slow',
                                                                              'deactivation (short pulse) t fast',
                                                                              'deactivation (short pulse) t mean']):
            self.plot_violin(param, labe)

    def plot_dea_long(self):

        for param, labe in zip(['dal_A1', 'dal_A2', 'dal_t1', 'dal_t2', 'dal_m'], ['deactivation (long pulse) A slow',
                                                                                   'deactivation (long pulse) A fast',
                                                                                   'deactivation (long pulse) t slow',
                                                                                   'deactivation (long pulse) t fast',
                                                                                   'deactivation (long pulse) t mean']):
            self.plot_violin(param, labe)


sns.set_style()
sns.set_context("talk")

metaWT = MetaWT()

metaWT.exp_lc_test()

metaWT.plot_rt()
metaWT.plot_fr()
metaWT.plot_ds()
metaWT.plot_ds1()
metaWT.plot_dea_long()
metaWT.plot_dea_short()
