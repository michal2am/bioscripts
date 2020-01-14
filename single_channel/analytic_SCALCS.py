import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import itertools as itr

from scalcs import scalcslib as scl
from scalcs import scplotlib as scpl
from scalcs import mechanism


class Model:

    def __init__(self, rates, resolution):
        """
        creates single model using SCALCS also handles analysis and printing/ploting for one model
        TODO: builds ony RFO SCALCS models, should be universal
        :param rates: list of rates [delta, gamma, beta, alpha]
        :param resolution: recording resolution for missed events correction
        """

        self.delta, self.gamma, self.beta, self.alpha = rates
        self.res = resolution
        self.mec = self.build_model()

        tres = self.res
        conc = 100e-9
        self.mec.set_eff('c', conc)
        results = scl.printout_distributions(self.mec, tres).split()

        to_frame = {'delta': [self.delta], 'gamma': [self.gamma], 'beta': [self.beta], 'alpha': [self.alpha], 'resolution': [self.res],
                    't1': results[16], 'a1': results[17], 't2': results[21], 'a2': results[22],
                    't1_a': results[50], 'a1_a': results[75], 't2_a': results[46], 'a2_a': results[73]
                    }

        self.frame = pd.DataFrame(to_frame)

    def build_model(self):
        """
        build RFO model
        TODO: model topology shall be specified
        :return: SCALCS mechanism model
        """

        mectitle = 'R-F-O'
        ratetitle = 'some numbers'
        O = mechanism.State('A', 'F*', 50e-12)
        F = mechanism.State('B', 'F', 0.0)
        R = mechanism.State('C', 'R', 0.0)
        RateList = [
            mechanism.Rate(self.alpha, O, F, name='alpha', limits=[1e-15, 1e+7]),
            mechanism.Rate(self.beta, F, O, name='beta', limits=[1e-15, 1e+7]),
            mechanism.Rate(self.gamma, F, R, name='gamma', limits=[1e-15, 1e+7]),
            mechanism.Rate(self.delta, R, F, name='delta', limits=[1e-15, 1e+10])
        ]

        return mechanism.Mechanism(RateList, mtitle=mectitle, rtitle=ratetitle)

    def print_mechanism(self):
        """
        prints standard SCALCS mechanism printout
        """
        print(self.mec)

    def print_occupancy(self):
        """
        prints states occupancy
        """
        print(scl.printout_occupancies(self.mec, self.res))

    def print_shut_dis(self):
        """
        prints shut time distribution
        TODO: only shuts, openings are commented out in SCALCS due to single open state issue
        """
        print(scl.printout_distributions(self.mec, self.res))

    def plot_shut_dis(self):
        """
        plots shut time distribution
        """
        t, ipdf, epdf, apdf = scpl.shut_time_pdf(self.mec, self.res)
        plt.semilogx(t, ipdf, 'r--', t, epdf, 'b-', t, apdf, 'g-')
        plt.ylabel('fshut(t)')
        plt.xlabel('Shut time, ms')
        plt.title('The shut time pdf')
        print('RED- ideal distribution\nGREEN- HJC distribution (corrected for missed events)')


class ModelRateGenerator:

    def __init__(self):

        self.deltas = self.rate_range_auto()
        self.gammas = self.rate_range_auto()
        self.betas = self.rate_range_auto()
        self.alphas = self.rate_range_auto()

        self.rate_sets = self.rate_combinations([self.deltas, self.gammas, self.betas, self.alphas])
        self.models = self.generate()

    def rate_range(self, start_rate):
        # step = int(start_rate * 0.did0)
        # return list(range(start_rate - 9 * step, start_rate + 10 * step, step))
        # return list(np.li100, 2500, 20))
        pass

    def rate_range_auto(self):
        return [int(rate) for rate in list(np.logspace(2.5, 4, num=15))]

    def rate_combinations(self, rates):
        return list(itr.product(*rates))

    def generate(self):

        models = []

        for rate_set in self.rate_sets:

            cfo_mechanism = Model(rate_set, 50e-6)
            models.append(cfo_mechanism.frame)

        all_models = pd.concat(models, axis=0)
        all_models.reset_index(inplace=True, drop=True)
        all_models.to_csv('all_test_log_SCALCS.csv')

        return all_models.melt(id_vars=['delta', 'gamma', 'beta', 'alpha']).copy()

class ModelPlots:

    def __init__(self, models):

        sns.set_style("white")
        sns.set_context("talk")

        self.models = models

    def correlation_plot(self):

        correlations_pearson = self.models.corr()
        annots_pearson = correlations_pearson.iloc[5:, 0:4].round(decimals=1)
        annots_pearson.replace(to_replace=[0.0, -0.0], value='', inplace=True)

        correlations_kendall = self.models.corr(method='kendall')
        annots_kendall = correlations_kendall.iloc[5:, 0:4].round(decimals=1)
        annots_kendall.replace(to_replace=[0.0, -0.0], value='', inplace=True)

        correlations_spearman = self.models.corr(method='spearman')
        annots_spearman = correlations_spearman.iloc[5:, 0:4].round(decimals=1)
        annots_spearman.replace(to_replace=[0.0, -0.0], value='', inplace=True)

        fig, axs = plt.subplots(ncols=3, nrows=1, figsize=[15, 5], sharey=True)
        fig.set_tight_layout(True)
        axs[0].set_title('pearson')
        axs[1].set_title('kendall')
        axs[2].set_title('spearman')
        sns.heatmap(correlations_pearson.iloc[5:, 0:4], yticklabels=1, cmap="RdBu_r", center=0, annot=annots_pearson, fmt='', ax=axs[0])
        sns.heatmap(correlations_kendall.iloc[5:, 0:4], yticklabels=1, cmap="RdBu_r", center=0, annot=annots_kendall, fmt='', ax=axs[1])
        sns.heatmap(correlations_spearman.iloc[5:, 0:4], yticklabels=1, cmap="RdBu_r", center=0, annot=annots_spearman, fmt='', ax=axs[2])

        fig.savefig('correlation')

    def one_rate_plots(self, property):

        rates = ['delta', 'gamma', 'beta', 'alpha']

        for rate in rates:

            fig, axs = plt.subplots(ncols=1, nrows=4, figsize=[10, 10], sharex=True)
            fig.set_tight_layout(True)
            hue_rates = rates.copy()
            hue_rates.remove(rate)

            sns.boxplot(x=rate, y=property, data=self.models, ax=axs[0], color='grey')

            for hue_rate, ax in zip(hue_rates, axs[1:]):
                sns.stripplot(x=rate, y=property, data=self.models, ax=ax, marker='.', color='red', size=3,
                                    hue=hue_rate, jitter=True, dodge=True)
                ax.get_legend().remove()
                ax.set_ylabel(property + ' | ' + hue_rate)
                ax.set_xticklabels(ax.get_xticklabels(), rotation='vertical')
            sns.despine(fig=fig)

            fig.savefig(property + '_' + rate)



new = Model([4300, 5000, 12500, 800], 50e-6)
new.plot_shut_dis()

#generate = ModelRateGenerator()
models = pd.read_csv('all_test_log_SCALCS.csv', index_col=0)

ploter = ModelPlots(models)

ploter.correlation_plot()
ploter.one_rate_plots('a2_a')

