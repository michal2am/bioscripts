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
        :param models: ModelBuilder object
        :param t0: simulation start time
        :param te: simulation end time
        :return:
        """
        self.models = models.models
        self.agonist_concentrations = models.agonist_concentrations
        self.stimuli = models.stimuli
        self.t0, self.te = t0, te

        # model parameters
        self.states = models.states
        self.states_names = models.states_names
        self.states_categories = models.states_categories
        self.states_belongs = models.states_belongs
        self.states_number = models.states_number
        self.states_ini_concentrations = models.states_ini_concentrations

        # integration results
        log.info('### Integration begins')
        self.results_dynamic = self.integrate_model(False)
        self.results_steady = self.integrate_model(True)

        # parsed results
        log.info('### Integration results parsing begins')
        self.tp_bystate, self.tp_bycategory = self.trajectories()
        self.ta_stimuli = self.trajectories_stimuli()
        self.steady_bystate, self.steady_bycategory = self.steady_occupancies()

        # global plot settings
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
            data.reset_index(inplace=True)
            data.loc[:, xlabel] = data.loc[:, xlabel].apply(np.log10)
            data.set_index(xlabel, inplace=True, drop=True)

            # prepare extended index to interpolate
            log_xknown = np.log10(xknown)
            log_xunknown = np.linspace(log_xknown[0], log_xknown[-1], num=100)
            log_x = np.concatenate((log_xknown[1:-1], log_xunknown), axis=0)
            log_x = np.sort(log_x)

            # insert new index and interpolate
            data = data.reindex(log_x)
            data = data.interpolate(method='akima')

            # invert index back
            data.reset_index(inplace=True)
            data.loc[:, xlabel] = data.loc[:, xlabel].rpow(np.ones(len(data[xlabel]))*10)
            data.set_index(xlabel, inplace=True, drop=True)

            return data

    def integrate_model(self, equi):
        """
        runs ode solver for all models
        :param equi:
        :return:
        """

        return [SolverOde(model.trmn, self.states_ini_concentrations, self.states_names, self.t0, self.te, equi) for model in self.models]

    def trajectories(self):
        """
        parse results from ode solver
        :return: list (concentrations) of dataframes (time + states/categories)
        """
        log.info('Parsing trajectories')
        tp_bystate, tp_bycategory = [], []
        for idx, result in enumerate(self.results_dynamic):
            tp_bystate.append(result.get_results())
            tp_bycategory.append(pd.DataFrame({category: tp_bystate[idx][self.states_belongs[category]].sum(axis=1) for category in self.states_belongs}))
        return [tp_bystate, tp_bycategory]

    def trajectories_stimuli(self):
        """
        generate stimulus trajectory
        :return: list (concentrations) of dataframes (time + stimuli)
        """
        log.info('Parsing stimuli')
        ta_stimuli = []
        for stimulus in self.stimuli:
            time = np.linspace(self.t0 - 2, self.te, int(1e5))
            df = pd.DataFrame({'time': time, 'stimuli': [stimulus(ti) for ti in time]})
            ta_stimuli.append(df.set_index('time', drop=True))

        return ta_stimuli

    def steady_occupancies(self):
        """
        parse steady state occupancies
        :return: list (concentrations) of dataframes (concentration + state/category)
        """
        log.info('Parsing steady state occupancies')
        steady_bystate, steady_bycategory = [], []
        for concentration, result in zip(self.agonist_concentrations, self.results_steady):
            steady = result.get_results().iloc[-1:]
            steady.index = [concentration]              # explicit to avoid nan substitution
            steady.index.name = 'concentration'
            steady_bystate.append(steady)
            steady_bycategory.append(pd.DataFrame({category: steady.loc[:, self.states_belongs[category]].sum(axis=1) for category in self.states_belongs}))

        return steady_bystate, steady_bycategory

    def plot_dynamic_response(self, concentration_index):
        """
        plots complete panel for selected concentration index
        :param concentration_index:
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
        self.tp_bystate[concentration_index].plot(ax=ax1, legend=True)
        self.tp_bycategory[concentration_index].plot(y='open',  ax=ax1, color='k', legend=True)

        ax1.legend(loc='upper center', ncol=5, borderaxespad=0.)
        ax1.set_xlim([self.t0 - 0.5, self.te + 0.5])
        ax1.set_ylim([0, 1.00])
        ax1.set_ylabel('state probability')

        # ax2 stimulus
        self.ta_stimuli[concentration_index].plot(ax=ax2, drawstyle='steps')

        ax2.set_xlim([self.t0 - 5, self.te + 2])
        ax2.set_ylim([0, self.agonist_concentrations[concentration_index]])
        ax2.set_yticks([0, self.agonist_concentrations[concentration_index]])
        ax2.set_ylabel('stimulus [mM]')
        ax2.set_xlabel('time [ms]')

        # ax3 equilibrium
        sns.barplot(data=self.steady_bystate[concentration_index], ax=ax3)

        log.info('###Equilibrium by states:')
        log.info(self.steady_bystate[concentration_index])

        ax3.set_ylabel('equilibrium occupancies')
        ax3.set_xlabel('state')
        ax3.set(yscale='log')

        # ax4 equilibrium categorical
        categories = self.steady_bycategory[concentration_index].transpose()
        categories.columns = ['cumulative equilibrium occupancy']
        categories.plot.pie(y='cumulative equilibrium occupancy', legend=True, labels=None, ax=ax4)

        log.info('###Equilibrium by categories:')
        log.info(self.steady_bycategory[concentration_index])

        ax4.legend(loc='lower center', nrow=1)  # throws no label warning, but it is ok
        ax4.axis('equal')

        sns.despine()

    def plot_steady_dose_response(self):
        """
        plots complete dose response curve
        :return:
        """

        fig_ode_dr = plt.figure()
        fig_ode_dr.canvas.set_window_title('ODE dose response results')
        gs = grd.GridSpec(2, 2, width_ratios=[1, 1], hspace=.1, wspace=0.)
        ax2 = plt.subplot(gs[1, 0:])
        ax1 = plt.subplot(gs[0, 0:], sharex=ax2)
        fig_ode_dr.add_subplot(ax1, ax2)

        # ax1 dr by state
        dr_bystate = pd.concat(self.steady_bystate)
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
        dr_bycategory = pd.concat(self.steady_bycategory)
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

