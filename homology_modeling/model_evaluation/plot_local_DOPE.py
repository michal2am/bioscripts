import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class PlotScopeEvaluateLocal:

    def __init__(self):

        self.singles = pd.read_csv("single.csv", index_col="residue number")
        self.sums = pd.read_csv("sumy.csv", index_col="chain")

    def plot_profile(self):
        # not normalized, each value has separate row
        all_singles = self.singles.drop(["mod_avg", "subunit"], axis=1)
        all_singles['template'], all_singles['type'], all_singles['subunit'] = all_singles['full_name'].str.split('#', 2).str
        all_singles = all_singles.drop(["full_name"], axis=1)
        all_singles = all_singles.melt(id_vars=["residue", "subunit", "position", "template_ene", 'template', 'type'])
        all_singles.rename(columns={'value': "absolute DOPE", 'variable': 'model'}, inplace=True)

        g = sns.FacetGrid(all_singles, col='template', row='type', col_order=['6i53', '6huk', '6hup', '6huo', '6huj', '6hug'],
                          row_order=['WT', 'GLY', 'LYS'], aspect=1, height=1.7125)
        g.map(sns.lineplot, "position", "absolute DOPE", "subunit", ci='sd')
        g.map(sns.lineplot, "position", "template_ene", "subunit", palette=sns.color_palette("pastel", 2))
        g.set_titles("{col_name}")

        #g.set(xlabel="residue", ylabel="absolute DOPE")
        g.set_axis_labels("residue", "absolute DOPE")
        plt.tight_layout()

        #g.fig.get_children()[-1].set_bbox_to_anchor((1.008, 1.15, 0, 0))

        #g.despine(trim=True)
        g.savefig("profile_DOPE.png", dpi=300)

    def plot_sum(self):

        # at this stage, each column is for one model, 4 rows for two chains 1 and 2 and DOPE absolute and relative if one receptor type
        n_sums = self.sums.iloc[:, :].div(self.sums.template_ene, axis=0)
        all_sums = pd.concat([self.sums, n_sums], axis=0)
        all_sums = all_sums.drop(["template_ene", "mod_avg"], axis=1)

        # as index had mixed category data, now those are as additional columns
        all_sums['template'], all_sums['type'], all_sums['subunit'] = all_sums.index.str.split('#', 2).str
        all_sums.reset_index(drop=True)
        # tidy up categories for relative values
        rows = int(all_sums.shape[0] / 2)
        all_sums['DOPE'] = ['absolute'] * rows + ['relative'] * rows

        # melting to have each measurement in a row, it gives 400 rows, 4 energies for 100 models if one receptor type
        all_sums = all_sums.melt(id_vars=['subunit', 'template', 'DOPE', 'type'])
        all_sums.rename(columns={'variable': 'model'}, inplace=True)
        all_sums["model"] = all_sums["model"].str.slice_replace(0, 6, '')

        # pivoting to have all categories as index and two columns of DOPE abs and rel score
        # this give most readable format
        all_sums = all_sums.pivot_table(index=['template', 'type', 'subunit', 'model'], columns='DOPE', values='value')
        print(all_sums)

        # reseting to put categories from index into columns, just for easier plotting with seaborn
        all_sums.reset_index(inplace=True)
        sorted = all_sums.sort_values(['type', 'subunit', 'relative', ], ascending=[True, True, False])
        sorted.to_csv('parsed_sums.csv')

        g1 = sns.catplot(x="type", y="absolute", hue="subunit", col="template",
                         kind="violin", inner="quart", split=True, data=all_sums, aspect=1, height=2.0,
                         col_order=['6i53', '6huk', '6hup', '6huo', '6huj', '6hug'], order=['WT', 'GLY', 'LYS'],
                         legend=False, col_wrap=4, sharex=False)
        g1.set_axis_labels("", "absolute DOPE")
        g1.set_titles("{col_name}")
        g1.set(ylim=(-0.61, -0.49))
        g1.despine(trim=True)
        plt.tight_layout()

        g1.savefig("violin_absoluteDOPE.png", dpi=150)

        g2 = sns.catplot(x="type", y="relative", hue="subunit", col="template",
                         kind="violin", inner="quart", split=True, data=all_sums, aspect=1, height=2.0,
                         col_order=['6i53', '6huk', '6hup', '6huo', '6huj', '6hug'], order=['WT', 'GLY', 'LYS'],
                         legend=False, col_wrap=4, sharex=False)
        g2.set_axis_labels("", "relative DOPE")
        g2.set_titles("{col_name}")
        g2.set(ylim=(0.8, 1.0))
        g2.despine(trim=True)
        plt.tight_layout()

        g2.savefig("violin_relativeDOPE.png", dpi=150)
'''
        g3 = sns.catplot(x="model", y="relative", hue="subunit", row='type', sharex='row', aspect=3,
                         kind="bar", data=all_sums)
        g3.set(ylim=(0.85, 1))
        g3.set_xticklabels(rotation=90)
        g3.set_xticklabels(fontsize=10)
        g3.fig.get_children()[-1].set_bbox_to_anchor((0.9, 1.05, 0, 0))

        g3.savefig("bar_relativeDOPE.png", dpi=300)
'''

if __name__ == '__main__':

    sns.set_style()
    sns.set_context("paper")

    plot = PlotScopeEvaluateLocal()

    plot.plot_sum()
    plot.plot_profile()

    plt.tight_layout()
    plt.show()
