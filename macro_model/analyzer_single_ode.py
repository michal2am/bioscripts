from solver_single_ode import SolverOde
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import seaborn as sns
import pandas as pd
import logging as log


class AnalyzerODE:

    def __init__(self, kinetic, t0, te):
        """
        :param kinetic: ModelBuilder object
        :param t0: simulation start time
        :param te: simulation end time
        :return:
        """
        self.kinetic = kinetic
        self.model = kinetic.models
        self.modeli_stimuli_idx = 0

        self.stimuli = kinetic.stimuli
        self.t0, self.te = t0, te

        # integration results
        self.results_dynamic = self.integrate_model()

        # parsed results
        self.tp_bystate, self.tp_bycategory = self.trajectories()
        self.ta_stimuli = self.trajectories_stimuli()

        # global plot settings
        sns.set_style("ticks", {'legend.frameon': True})
        sns.set_context("poster")
        sns.set_palette('Paired', kinetic.states_number)
        self.colors = sns.color_palette('Paired', kinetic.states_number)

    def integrate_model(self):
        """
        runs ode solver for all models
        :return:
        """
        log.info('### Integration begins ###')
        return SolverOde(self.model[self.modeli_stimuli_idx].trmn, self.kinetic.states_ini_concentrations, self.kinetic.states_names, self.t0, self.te)

    def trajectories(self):
        """
        parse results from ode solver
        :return: list (concentrations) of dataframes (time + states/categories)
        """
        log.info('### Parsing trajectories ###')
        return [self.results_dynamic.get_results(),
                pd.DataFrame({category: self.results_dynamic.get_results()[self.kinetic.states_belongs[category]].sum(axis=1) for category in self.kinetic.states_belongs})]

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

    def plot_dynamic_response(self):
        """
        plots complete panel for selected concentration index
        :return:
        """

        # set grid
        fig_ode = plt.figure()
        fig_ode.canvas.set_window_title('ODE results')
        gs = grd.GridSpec(6, 2, width_ratios=[3, 2], hspace=2, wspace=0.50)
        ax2 = plt.subplot(gs[5, 0])
        ax1 = plt.subplot(gs[0:5, 0], sharex=ax2)
        fig_ode.add_subplot(ax1, ax2)

        # ax1 trajectory

        # reverse values and set params
        self.tp_bystate.apply(lambda x: x*-1).plot(ax=ax1, legend=True)
        self.tp_bycategory.apply(lambda x: x*-1).plot(y='open', ax=ax1, linestyle=':', color='grey', linewidth=5, legend=False)

        # additional editing
        ax1.legend(bbox_to_anchor=(1.05, 1), ncol=3, borderaxespad=0.)
        ax1.set_xlim([self.t0 - 0.5, self.te + 0.5])
        ax1.set_ylim([-100.00, 0.1])
        ax1.set_ylabel('state probability')

        # ax2 stimulus

        # get stimuli values
        self.ta_stimuli[self.modeli_stimuli_idx].plot(ax=ax2, drawstyle='steps')

        # additional editing
        ax2.legend(bbox_to_anchor=(1.05, 1))
        ax2.set_xlim([self.t0 - 5, self.te + 2])
        ax2.set_ylim([0, self.kinetic.agonist_concentrations[self.modeli_stimuli_idx]])
        ax2.set_yticks([0, self.kinetic.agonist_concentrations[self.modeli_stimuli_idx]])
        ax2.set_ylabel('stimulus [mM]')
        ax2.set_xlabel('time [ms]')

        sns.despine()
