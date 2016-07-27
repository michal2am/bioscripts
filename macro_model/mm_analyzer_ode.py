from mm_solver_ode import SolverOde
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import matplotlib.colors as cls
import seaborn as sns
import pandas as pd
import logging as log


class AnalyzerODE:

    def __init__(self, models, t0, te):
        """

        :param models:
        :param t0:
        :param te:
        :return:
        """
        self.models = models.models
        self.agonist_concentrations = models.agonist_concentrations
        self.stimuli = models.stimuli
        self.t0, self.te = t0, te

        self.states = models.states
        self.states_names = models.states_names
        self.states_categories = models.states_categories
        self.states_belongs = models.states_belongs
        self.states_number = models.states_number
        self.states_ini_concentrations = models.states_ini_concentrations

        self.results_dynamic = self.integrate_model(False)
        self.results_steady = self.integrate_model(True)

        sns.set_style("ticks", {'legend.frameon': True})
        sns.set_context("poster")
        sns.set_palette('Paired', self.states_number)
        self.colors = sns.color_palette('Paired', self.states_number)

    @staticmethod
    def fit_logarithmic(data):
            """
            fits numeric interpolation with logarithmic x-values conversion
            :param data: pandas dataframe, x-values as index
            :return: pandas dataframe, interpolated x-values as index
            """

            xlabel = data.index.name
            xknown = data.index

            # replace index with index logarithm
            data = data.reset_index()
            data[xlabel] = data[xlabel].apply(np.log10)
            data = data.set_index(xlabel, drop=True)

            # prepare extended index to interpolate
            log_xknown = np.log10(xknown)
            log_xunknown = np.linspace(log_xknown[0], log_xknown[-1], num=100)
            log_x = np.concatenate((log_xknown[1:-1], log_xunknown), axis=0)
            log_x = np.sort(log_x)

            # insert new index and interpolate
            data = data.reindex(log_x)
            data = data.interpolate(method='akima')

            # invert index back
            data = data.reset_index()
            data[xlabel] = data[xlabel].rpow(np.ones(len(data[xlabel]))*10)
            data = data.set_index(xlabel, drop=True)

            return data

    def integrate_model(self, equi):
        """

        :param equi:
        :return:
        """
        return [SolverOde(model.trmn, self.states_ini_concentrations, self.states_names, self.t0, self.te, equi) for model in self.models]

    def dynamic_response(self):
        """

        :return:
        """

        fig_ode = plt.figure()
        fig_ode.canvas.set_window_title('ODE results')
        gs = grd.GridSpec(6, 2, width_ratios=[3, 2], hspace=2, wspace=0.20)
        ax2 = plt.subplot(gs[5, 0])
        ax1 = plt.subplot(gs[0:5, 0], sharex=ax2)
        ax4 = plt.subplot(gs[3:, 1])
        ax3 = plt.subplot(gs[0:3, 1])
        fig_ode.add_subplot(ax1, ax2, ax3, ax4)

        # ax1 trajectory
        tp_bystate = self.results_dynamic[0].get_results()
        tp_bystate.plot(x='time', ax=ax1, legend=True)

        tp_bycategory = pd.DataFrame({category: tp_bystate.drop('time', 1)[self.states_belongs[category]].sum(axis=1) for category in self.states_belongs})
        tp_bycategory = pd.concat([tp_bystate['time'], tp_bycategory], axis=1)
        tp_bycategory.plot(x='time',  y='open',  ax=ax1, color='k', legend=True)

        ax1.legend(loc='upper center', ncol=5, borderaxespad=0.)
        ax1.set_xlim([self.t0 - 0.5, self.te + 0.5])
        ax1.set_ylim([0, 1.00])
        ax1.set_ylabel('state probability')

        # ax2 stimulus
        time = np.linspace(self.t0 - 2, self.te, 1e5)
        stim = [self.stimuli[0](ti) for ti in time]
        ax2.step(time, stim, 'k-')

        ax2.set_xlim([self.t0 - 5, self.te + 2])
        ax2.set_ylim([0, max(stim)])
        ax2.set_yticks([0, max(stim)])
        ax2.set_ylabel('stimulus [mM]')
        ax2.set_xlabel('time [ms]')

        # ax3 equilibrium
        steady_bystate = self.results_steady[0].get_results()[-1:].drop('time', 1)
        log.info('###Equilibrium by states:')
        log.info(steady_bystate)
        sns.barplot(data=steady_bystate, ax=ax3)
        ax3.set_ylabel('equilibrium occupancies')
        ax3.set_xlabel('state')
        ax3.set(yscale='log')

        # ax4 equilibrium categorical
        steady_bycategory = pd.DataFrame({category: steady_bystate[self.states_belongs[category]].sum(axis=1) for category in self.states_belongs})
        log.info('###Equilibrium by categories:')
        log.info(steady_bycategory)
        categories = steady_bycategory.transpose()
        categories.columns = ['cumulative equilibrium occupancy']
        categories.plot.pie(y='cumulative equilibrium occupancy', legend=True, labels=None, ax=ax4)
        ax4.legend(loc='lower center', nrow=1)
        ax4.axis('equal')

        sns.despine()

    def steady_dose_response(self):
        """

        :return:
        """

        fig_ode_dr = plt.figure()
        fig_ode_dr.canvas.set_window_title('ODE dose response results')
        gs = grd.GridSpec(2, 2, width_ratios=[1, 1], hspace=.1, wspace=0.)
        ax2 = plt.subplot(gs[1, 0:])
        ax1 = plt.subplot(gs[0, 0:], sharex=ax2)
        fig_ode_dr.add_subplot(ax1, ax2)

        # ax1 dr by state
        concs = pd.DataFrame({'concentration': self.agonist_concentrations})
        dr_bystate = pd.DataFrame()

        for result in self.results_steady:
            dr_bystate = dr_bystate.append(result.get_results()[-1:].drop('time', 1))

        dr_bystate = dr_bystate.reset_index(drop=True)
        dr_bystate = pd.concat([concs, dr_bystate], axis=1)
        dr_bystate = dr_bystate.set_index('concentration', drop=True)
        dr_bystate_fit = dr_bystate.copy()
        dr_bystate_fit = self.fit_logarithmic(dr_bystate_fit)
        dr_bystate.plot(ax=ax1, marker='o', linestyle='None', legend=False)
        dr_bystate_fit.plot(ax=ax1, legend=False)

        log.info('###Equilibrium by states:')
        log.info(dr_bystate)

        ax1.legend(labels=self.states_names, numpoints=1, loc='center right', ncol=2, borderaxespad=0.)
        ax1.set(xscale='log')
        ax1.set_ylabel("equilibrium occupancy")

        # ax2 dr bystate
        dr_bycategory = pd.DataFrame({category: dr_bystate[self.states_belongs[category]].sum(axis=1) for category in self.states_belongs})
        dr_bycategory_fit = dr_bycategory.copy()
        dr_bycategory_fit = self.fit_logarithmic(dr_bycategory_fit)
        dr_bycategory.plot(ax=ax2, style='o', colormap=cls.ListedColormap(self.colors[0:6]))
        dr_bycategory_fit.plot(ax=ax2, legend=None, colormap=cls.ListedColormap(self.colors[0:6]))

        log.info('###Equilibrium by categories:')
        log.info(dr_bycategory)

        ax2.legend(labels=self.states_categories, numpoints=1, loc='center right', ncol=2, borderaxespad=0.)
        ax2.set(xscale='log')
        ax2.set_xlabel("stimulus [mM]")
        ax2.set_ylabel("equilibrium occupancy")
        ax2.set_xlim(1e-4, 1e3)

        sns.despine()



