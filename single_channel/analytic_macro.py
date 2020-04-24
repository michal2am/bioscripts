# script for some trend-analytic single channel simulations
# if "generate" is commented out simulations results are read from csv

# TODO: a lot of stuff, thus is just temporary (see detailed TODOS)

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import itertools as itr

from scalcs import scalcslib as scl
from scalcs import scplotlib as scpl
from scalcs import mechanism


class ModelRFO:

    def __init__(self, rates):
        """

        :param rates:
        """

        self.mechanism = rates
        print('Building new model')

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

    def __init__(self):
        """
        creates mutliple possible models
        TODO: works only for CFO model
        TODO: works only for logarithmic parameter sets
        """

        self.deltas = self.rate_range_auto()
        self.gammas = self.rate_range_auto()
        self.betas = self.rate_range_auto()
        self.alphas = self.rate_range_auto()

        self.rate_sets = self.rate_combinations([self.deltas, self.gammas, self.betas, self.alphas])
        self.models = self.generate()

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
        return [int(rate) for rate in list(np.logspace(2.7, 4.31, num=25))]

    def rate_combinations(self, rates):
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

            cfo_mechanism = ModelAnalysis('RFO', rate_set, 50e-6)
            generated_models.append(cfo_mechanism.frame)

        all_models = pd.concat(generated_models, axis=0)
        all_models.reset_index(inplace=True, drop=True)
        all_models.to_csv('all_test_log_SCALCS_REDONE_25bis.csv')

        return all_models.melt(id_vars=['delta', 'gamma', 'beta', 'alpha']).copy()


generate = ModelMultiGenerator()
#models = pd.read_csv('all_test_log_SCALCS_REDONE_25bis.csv', index_col=0)

