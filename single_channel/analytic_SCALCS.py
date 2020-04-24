# script for some trend-analytic single channel simulations
# if "generate" is commented out simulations results are read from csv

# TODO: a lot of stuff, thus is just temporary (see detailed TODOS)

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import itertools as itr
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
import argparse


from scalcs import scalcslib as scl
from scalcs import scplotlib as scpl
from scalcs import mechanism


class ModelRFO:

    def __init__(self, rates):
        """

        :param rates:
        """

        self.mechanism = rates

    @property
    def mechanism(self):
        return self._mechanism

    @mechanism.setter
    def mechanism(self, rates):

        o = mechanism.State('A', 'F*', 50e-12)
        f = mechanism.State('B', 'F', 0.0)
        r = mechanism.State('C', 'R', 0.0)

        rate_list = [
            mechanism.Rate(rates[3], o, f, name='alpha', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates[2], f, o, name='beta', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates[1], f, r, name='gamma', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates[0], r, f, name='delta', limits=[1e-15, 1e+10])
        ]

        complete_mechanism = mechanism.Mechanism(rate_list, mtitle='CFO', rtitle='CFO_rates')
        complete_mechanism.set_eff('c', 100e-9)

        self._mechanism = complete_mechanism


class ModelAnalysis:

    def __init__(self, topology, rates, resolution):
        """
        creates single model using SCALCS also handles analysis and printing/ploting for one model
        :param topology:
        :param rates: list of rates
        :param resolution: recording resolution for missed events correction
        """

        self.resolution = resolution
        self.topology = topology
        self.mechanism = rates

        # TODO: clean up parameter calculation and data frame build, should work for any model
        # TODO: dirty ideal open only for RFO
        parameters = {rate.name: rate.rateconstants[0] for rate in self.mechanism.Rates}
        parameters.update({'resolution': self.resolution})
        shuts_dist = scl.printout_distributions(self.mechanism, self.resolution).split()
        parameters.update({'t1': shuts_dist[16], 'a1': shuts_dist[17], 't2': shuts_dist[21], 'a2': shuts_dist[22],
                        't1_a': shuts_dist[50], 'a1_a': shuts_dist[75], 't2_a': shuts_dist[46], 'a2_a': shuts_dist[73],
                        'op_t': 1/self.mechanism.Rates[0].rateconstants[0]*1000})

        self.frame = pd.DataFrame([parameters])

    @property
    def topology(self):
        return self._topology

    @topology.setter
    def topology(self, topology):
        self._topology = topology

    @property
    def resolution(self):
        return self._res

    @resolution.setter
    def resolution(self, resolution):
        self._res = resolution

    @property
    def mechanism(self):
        return self._mechanism

    @mechanism.setter
    def mechanism(self, rates):
        if self.topology == 'RFO':
            self._mechanism = ModelRFO(rates).mechanism

    def print_mechanism(self):
        """
        prints standard SCALCS mechanism printout
        """
        print(self.mechanism)

    def print_occupancy(self):
        """
        prints states occupancy
        """
        print(scl.printout_occupancies(self.mechanism, self.resolution))

    def print_shut_dis(self):
        """
        prints shut time distribution
        TODO: only shuts, openings are commented out in SCALCS due to single open state issue
        """
        print(scl.printout_distributions(self.mechanism, self.resolution))

    def plot_shut_dis(self):
        """
        plots shut time distribution
        """
        t, ipdf, epdf, apdf = scpl.shut_time_pdf(self.mechanism, self.resolution)
        plt.semilogx(t, ipdf, 'r--', t, epdf, 'b-', t, apdf, 'g-')
        plt.ylabel('fshut(t)')
        plt.xlabel('Shut time, ms')
        plt.title('The shut time pdf')


class ModelMultiGenerator:

    def __init__(self, resolution, out_table_file):
        """
        creates mutliple possible models
        TODO: works only for CFO model
        TODO: works only for logarithmic parameter sets
        """

        self.resolution = resolution
        self.output = out_table_file

        # TODO: all below is a mess
        self.deltas = self.rate_range_auto()
        self.gammas = self.rate_range_auto()
        self.betas = self.rate_range_auto()
        self.alphas = self.rate_range_auto()

        self.rate_sets = self.rate_combinations([self.deltas, self.gammas, self.betas, self.alphas])
        self.models = self.generate()

    @property
    def resolution(self):
        return self._resolution

    @resolution.setter
    def resolution(self, resolution):
        self._resolution = resolution

    @property
    def output(self):
        return self._output

    @output.setter
    def output(self, output_table_file):
        self._output = output_table_file

    def rate_range(self, start_rate):
        """
        TODO: shuld be a part of custom rate range generation option
        """
        # step = int(start_rate * 0.did0)
        # return list(range(start_rate - 9 * step, start_rate + 10 * step, step))
        # return list(np.li100, 2500, 20))
        pass

    @staticmethod
    def rate_range_auto():
        """
        generates a logarithmic set of rate values for single rate
        :return: list of rates values
        """
        # return [int(rate) for rate in list(np.logspace(2.7, 4.31, num=25))]
        return [int(rate) for rate in list(np.logspace(2.4, 4.5, num=30))]

    @staticmethod
    def rate_combinations(rates):
        """
        generates all possible permutations of rates for complete model
        :param rates: list of list rates
        :return: list of lists of rate values sets
        """
        return list(itr.product(*rates))

    def generate(self):
        """
        generates multiple models based on rates sets
        :return: dataframe of models
        """

        generated_models = []

        for rate_set in self.rate_sets:

            resolution = self.resolution
            cfo_mechanism = ModelAnalysis('RFO', rate_set, resolution)
            generated_models.append(cfo_mechanism.frame)

        all_models = pd.concat(generated_models, axis=0)
        all_models.reset_index(inplace=True, drop=True)
        all_models.to_csv(self.output)

        return all_models.melt(id_vars=['delta', 'gamma', 'beta', 'alpha']).copy()


class ModelPlots:

    def __init__(self, models_data, rates_names):
        """

        :param models_data: pandas dataframe of model data, columns: rates / resolution / parameters
        :param rates_names: list of rates column names
        """

        sns.set_style("white")
        sns.set_context("talk")

        self.models = models_data
        self.rates = rates_names
        self.no_rates = rates_names

    @property
    def models(self):
        return self._models

    @models.setter
    def models(self, models_data):
        self._models = models_data

    @property
    def rates(self):
        return self._rates

    @rates.setter
    def rates(self, rates_names):
        self._rates = rates_names

    @property
    def no_rates(self):
        return self._no_rates

    @no_rates.setter
    def no_rates(self, rates_names):
        self._no_rates = len(rates_names)

    def correlation_plot(self):
        """
        plots correlation between rates and parameters
        """

        correlations_pearson = self.models.corr()
        annots_pearson = correlations_pearson.iloc[self.no_rates + 1:, 0:self.no_rates].round(decimals=1)
        annots_pearson.replace(to_replace=[0.0, -0.0], value='', inplace=True)

        correlations_kendall = self.models.corr(method='kendall')
        annots_kendall = correlations_kendall.iloc[self.no_rates + 1:, 0:self.no_rates].round(decimals=1)
        annots_kendall.replace(to_replace=[0.0, -0.0], value='', inplace=True)

        correlations_spearman = self.models.corr(method='spearman')
        annots_spearman = correlations_spearman.iloc[self.no_rates + 1:, 0:self.no_rates].round(decimals=1)
        annots_spearman.replace(to_replace=[0.0, -0.0], value='', inplace=True)

        fig, axs = plt.subplots(ncols=3, nrows=1, figsize=[15, 5], sharey=True)
        fig.set_tight_layout(True)
        axs[0].set_title('pearson')
        axs[1].set_title('kendall')
        axs[2].set_title('spearman')
        sns.heatmap(correlations_pearson.iloc[self.no_rates + 1:, 0:self.no_rates], yticklabels=1,
                    cmap="RdBu_r", center=0, annot=annots_pearson, fmt='', ax=axs[0])
        sns.heatmap(correlations_kendall.iloc[self.no_rates + 1:, 0:self.no_rates], yticklabels=1,
                    cmap="RdBu_r", center=0, annot=annots_kendall, fmt='', ax=axs[1])
        sns.heatmap(correlations_spearman.iloc[self.no_rates + 1:, 0:self.no_rates], yticklabels=1,
                    cmap="RdBu_r", center=0, annot=annots_spearman, fmt='', ax=axs[2])

        fig.savefig('correlation')

    def one_rate_plots(self, parameter):
        """
        plots selected parameter in function of rates,
        separate panel for each rate, first row cumulative boxplot, other rows stripps for single other rate dependency
        :param parameter: selected parameter
        """

        for rate in self.rates:

            fig, axs = plt.subplots(ncols=1, nrows=self.no_rates, figsize=[10, 10], sharex=True)
            fig.set_tight_layout(True)
            hue_rates = self.rates.copy()
            hue_rates.remove(rate)

            sns.boxplot(x=rate, y=parameter, data=self.models, ax=axs[0], color='grey')

            for hue_rate, ax in zip(hue_rates, axs[1:]):
                sns.stripplot(x=rate, y=parameter, data=self.models, ax=ax, marker='.', color='red', size=3,
                              hue=hue_rate, jitter=True, dodge=True)
                ax.get_legend().remove()
                ax.set_ylabel(parameter + ' | ' + hue_rate)
                ax.set_xticklabels(ax.get_xticklabels(), rotation='vertical')
            sns.despine(fig=fig)

            fig.savefig(parameter + '_' + rate)


class ModelParse:

    def __init__(self, fit, models_csvs, experiments_csv):

        if fit:

            self.models = models_csvs
            self.experiment = experiments_csv

        else:

            pass

    @property
    def models(self):
        return self._models

    @models.setter
    def models(self, model_files):
        models = []
        for model_file in model_files:
            model = pd.read_csv(model_file, index_col=0)
            models.append(model)
        models_df = pd.concat(models)
        self._models = models_df.drop(columns=['t1', 'a1', 't2', 'a2'])

    @property
    def experiment(self):
        return self._experiment

    @experiment.setter
    def experiment(self, experimental_table):
        experiment = pd.read_csv(experimental_table)
        self._experiment = experiment.drop(columns=['file', 'p1_s', 'p2_s', 't3_s', 'p3_s', 't4_s', 'p4_s', 'mean_s', 't1_o', 'p1_o', 't2_o' ,'p2_o'])

    def parse_and_fit(self):
        # TODO: refactor indexing, it's a mess
        """

        :return:
        """

        fitted = pd.DataFrame(columns=['type', 'project', 'equilibrium', 'forward'])

        for index, row in self.experiment.iterrows():

            experimental = row.to_frame().T.rename(columns={'t1_s': 't2_a', 't2_s': 't1_a', 'p1_n': 'a2_a', 'p2_n': 'a1_a', 'mean_o': 'op_t'})
            experimental = experimental.loc[:, ['type', 'project', 'resolution', 't1_a', 'a1_a', 't2_a', 'a2_a', 'op_t']]

            res_filter = self.models[self.models.loc[:, 'resolution'] == experimental.loc[:, 'resolution'].values[0]]

            up = 1.8
            bt = 0.2

            pre_filter = res_filter[(res_filter.t2_a < up*row.t1_s) & (res_filter.t2_a > bt*row.t1_s) &
                                 (res_filter.a2_a < up*row.p1_n) & (res_filter.a2_a > bt*row.p1_n) &
                                 (res_filter.t1_a < up*row.t2_s) & (res_filter.t1_a > bt*row.t2_s) &
                                 (res_filter.a1_a < up*row.p2_n) & (res_filter.a1_a > bt*row.p2_n) &
                                 (res_filter.op_t < 1.0*row.mean_o) & (res_filter.op_t > 0.2*row.mean_o)]

            relatives = {prop: pre_filter.loc[:, prop].apply(lambda x: abs((x - experimental.loc[index, prop])/x))
                         for prop in ['t1_a', 'a1_a', 't2_a', 'a2_a', 'op_t']}
            relatives = pd.DataFrame.from_dict(relatives)
            relatives.loc[:, 'diff_av'] = relatives.mean(numeric_only=True, axis=1)

            pre_filter = pd.concat([pre_filter, relatives.loc[:, 'diff_av']], axis=1)
            post_filter = pre_filter.sort_values(by=['diff_av']).iloc[0:1, :]

            post_filter.loc[:, 'forward'] = np.log(post_filter.loc[:, 'beta'] * post_filter.loc[:, 'delta'])
            post_filter.loc[:, 'equilibrium'] = np.log((post_filter.loc[:, 'beta'] * post_filter.loc[:, 'delta'])/
                                                       (post_filter.loc[:, 'alpha'] * post_filter.loc[:, 'gamma']))

            post_filter.reset_index(inplace=True, drop=True)
            experimental.reset_index(inplace=True, drop=True)

            print(experimental)
            print(post_filter)

            result = pd.concat([experimental[['type', 'project']], post_filter[['equilibrium', 'forward']]], axis=1)
            fitted = fitted.append(result, ignore_index=True)

        fitted.to_csv('rates.csv')

    def plot_phi(self):
        """

        :return: 
        """

        data = pd.read_csv('rates.csv', index_col=0)
        sns.set_style("white")
        sns.set_context("talk")

        phis = pd.DataFrame(columns=['project', 'phi_s', 'phi_c'])


        for project in data.project.unique():

            single_cell = data[data.loc[:, 'project'] == project]
            cumulative_cell = single_cell.groupby('type', sort=False).mean().reset_index()

            g1 = sns.relplot(x='equilibrium', y='forward', hue='type',  data=single_cell)
            g1.map(sns.regplot, x='equilibrium', y='forward', scatter=False, data=single_cell)
            plt.show()
            g1.savefig(project+'_single_REFER.png')

            g2 = sns.relplot(x='equilibrium', y='forward', hue='type', data=cumulative_cell)
            g2.map(sns.regplot, x='equilibrium', y='forward', scatter=False, data=cumulative_cell)
            plt.show()
            g2.savefig(project + '_cumulative_REFER.png')


            x_s = single_cell.loc[:, 'equilibrium']
            x_s = sm.add_constant(x_s)
            y_s = single_cell.loc[:, 'forward']
            model_s = sm.OLS(y_s, x_s).fit()

            x_c = cumulative_cell.loc[:, 'equilibrium']
            x_c = sm.add_constant(x_c)
            y_c = cumulative_cell.loc[:, 'forward']
            model_c = sm.OLS(y_c, x_c).fit()

            phis = phis.append({'project': project, 'phi_s': model_s.params.loc['equilibrium'],
                               'phi_c': model_c.params.loc['equilibrium']}, ignore_index=True)

        phis.to_csv('phis.csv')


parser = argparse.ArgumentParser()
parser.add_argument('--mode', help='foo help')
parser.add_argument('--generate_resolution', type=float, help='foo help')
parser.add_argument('--generate_output', help='foo help')
parser.add_argument('--experimental_input', help='foo help')


args = parser.parse_args()

#new = ModelAnalysis('RFO', [4300, 5000, 12500, 800], 50e-6)
#new.plot_shut_dis()
#new.print_mechanism()

if args.mode == 'generate':

    generate = ModelMultiGenerator(args.generate_resolution, args.generate_output)

if args.mode == 'plot_trend':

    pass
    # TODO: should accept file, not DataFrame
    # TODO: should work with multi resolution DataFrame(s) from file list

    # models = pd.read_csv('all_test_log_SCALCS_REDONE_25bis.csv', index_col=0)
    # ploter = ModelPlots(models, ['delta', 'gamma', 'beta', 'alpha'])
    # ploter.correlation_plot()
    # ploter.one_rate_plots('a2_a')

if args.mode == 'fit_plot':

    parser = ModelParse(True, ['all_test_log_SCALCS_30.csv', 'all_test_log_SCALCS_30_r30.csv',
                         'all_test_log_SCALCS_30_r40.csv', 'all_test_log_SCALCS_30_r60.csv',
                         'all_test_log_SCALCS_30_r70.csv', 'all_test_log_SCALCS_30_r80.csv'],
                        args.experimental_input)

    parser.parse_and_fit()
    parser.plot_phi()

if args.mode == 'fitted_plot':

    parser = ModelParse(False, None, None)
    parser.plot_phi()
