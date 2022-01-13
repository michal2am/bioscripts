import numpy as np
import pandas as pd
import plotly.express as px
from outliers import smirnov_grubbs as grubbs
from scikit_posthocs import outliers_iqr
import plotly.graph_objects as go
import statsmodels
import argparse
import re
from plotly.subplots import make_subplots


class REFERPloter:

    @staticmethod
    def forward_flip(cell):
        return np.log(cell['delta'])

    @staticmethod
    def equilibrium_flip(cell):
        return np.log(cell['delta'] / cell['gamma'])

    @staticmethod
    def forward_gate(cell):
        return np.log(cell['beta'])

    @staticmethod
    def equilibrium_gate(cell):
        return np.log(cell['beta'] / cell['alpha'])

    @staticmethod
    def energy_flip(cell):
        return -0.59 * np.log(cell['delta'] / cell['gamma'])

    @staticmethod
    def energy_gate(cell):
        return -0.59 * np.log(cell['beta'] / cell['alpha'])

    def outliers_grubbs(self):
        '''
        apparently something wrong, the test is oversensitive compared to graphpad
        TODO: make it work like graphpad ???
        :return:
        '''

        for mut in self.table_results.type.unique():
            for rate in ['alpha', 'beta', 'gamma', 'delta']:

                rates = list(
                    self.table_results.loc[self.table_results['type'] == mut][rate])  # list() just for index reset
                cells = list(self.table_results.loc[self.table_results['type'] == mut]['file'])
                print(rates)
                outliers_val = grubbs.two_sided_test_outliers(rates, alpha=.05)
                outliers_ind = grubbs.two_sided_test_indices(rates, alpha=.05)

                if outliers_val:
                    print(rates)
                    outliers_cel = [cells[ind] for ind in outliers_ind]
                    print('{} value outliers in {} type; cells: {}, values: {}'.format(rate, mut, outliers_cel,
                                                                                       outliers_val))
                    self.table_results.loc[self.table_results['file'].isin(outliers_cel), rate] = np.NAN

    def outliers_iqrsckit(self):

        for mut in self.table_results.type.unique():
            for rate in ['alpha', 'beta', 'gamma', 'delta']:

                rates = list(
                    self.table_results.loc[self.table_results['type'] == mut][rate])  # list() just for index reset
                cells = list(self.table_results.loc[self.table_results['type'] == mut]['file'])
                outliers_val = outliers_iqr(rates, ret='outliers')
                outliers_ind = outliers_iqr(rates, ret='outliers_indices')

                if outliers_val.size > 0:
                    outliers_cel = [cells[ind] for ind in outliers_ind]
                    print('{} value outliers in {} type; cells: {}, values: {} \n within: {}'.
                          format(rate, mut, outliers_cel, outliers_val, rates))
                    self.table_results.loc[self.table_results['file'].isin(outliers_cel), rate] = np.NAN

    def __init__(self, rates_table, project, config):

        self.config = config
        self.project = project

        self.table_results = pd.read_csv(rates_table).drop(
            columns=['Unnamed: 0', 'e', 'f1', 'f2', 'f3'])  # remove e and fx drops after fit script cleanup
        self.wt_results = pd.read_csv('hjcfit_rates_WT_pooled.csv').drop(columns=['Unnamed: 0', 'e', 'f1', 'f2', 'f3'])

        self.table_results = pd.concat([self.table_results, self.wt_results])     # pooled WT results
        self.table_results['project'] = self.project

        self.outliers_iqrsckit()

        self.cumulative_results = self.table_results.drop(columns=['file'])
        self.cumulative_results_m = self.cumulative_results.groupby(
            ['project', 'type', 'model'], as_index=False)[['alpha', 'beta', 'gamma', 'delta']].mean()

        self.cumulative_results_m['forward_flip'] = self.cumulative_results_m.apply(
            lambda cell: self.forward_flip(cell), axis=1)
        self.cumulative_results_m['equilibrium_flip'] = self.cumulative_results_m.apply(
            lambda cell: self.equilibrium_flip(cell), axis=1)
        self.cumulative_results_m['forward_gate'] = self.cumulative_results_m.apply(
            lambda cell: self.forward_gate(cell), axis=1)
        self.cumulative_results_m['equilibrium_gate'] = self.cumulative_results_m.apply(
            lambda cell: self.equilibrium_gate(cell), axis=1)
        self.cumulative_results_m['energy_flip'] = self.cumulative_results_m.apply(
            lambda cell: self.energy_flip(cell), axis=1)
        self.cumulative_results_m['energy_gate'] = self.cumulative_results_m.apply(
            lambda cell: self.energy_gate(cell), axis=1)

        # % errors
        self.cumulative_results_d = self.cumulative_results.groupby(['project', 'type', 'model'], as_index=False)[
            ['alpha', 'beta', 'gamma', 'delta']].apply(lambda x: np.std(x) / np.mean(x) * 100)

        self.table_results.sort_values('type', ascending=False, inplace=True)
        self.cumulative_results_m.sort_values('type', ascending=False, inplace=True)
        self.cumulative_results_d.sort_values('type', ascending=False, inplace=True)

        print(self.cumulative_results_m)
        # print(self.cumulative_results_d)

        self.cumulative_results_m.to_csv('hjcfit_rates_m_' + self.project + '.csv')
        self.cumulative_results_d.to_csv('hjcfit_rates_d_' + self.project + '.csv')

    def REFER_plot_Auerbach(self):

        plots = []
        for axes in zip(['phi_flip', 'phi_gate', 'efficiency'],
                        ['equilibrium_flip', 'equilibrium_gate', 'energy_gate'],
                        ['forward_flip', 'forward_gate', 'energy_flip']):
            plot = (px.scatter(self.cumulative_results_m, x=axes[1], y=axes[2],
                               title=axes[0],
                               color='type', template='presentation', width=400, height=400,
                               hover_name='type', hover_data=['alpha', 'beta', 'gamma', 'delta'],
                               color_discrete_sequence=px.colors.qualitative.Dark24,
                               ))
            plot.add_trace(
                px.scatter(self.cumulative_results_m, x=axes[1], y=axes[2],
                           trendline='ols',
                           color_discrete_sequence=px.colors.qualitative.Dark24,
                           ).data[1])

            plots.append(plot)

        with open(self.project + '_auerbach.html', 'w') as f:
            for plot in plots:
                f.write(plot.to_html(full_html=False, include_plotlyjs='cdn'))

    def REFER_plot_Neudecker(self):

        for step in ['f1', 'f2', 'f3']:
            project_allCells = (
                px.scatter(self.table_results[self.table_results['project'] == self.project], x='e', y=step,
                           title=self.project + 'single point - single cell',
                           color='type', template='presentation', width=8000, height=800,
                           marginal_x='rug', marginal_y='rug',
                           hover_name='file', hover_data=['alpha', 'beta', 'gamma', 'delta'],
                           color_discrete_sequence=px.colors.qualitative.Dark24,
                           ))
            project_allCells.add_trace(
                px.scatter(self.table_results[self.table_results['project'] == self.project], x='e', y=step,
                           trendline='ols',
                           color_discrete_sequence=px.colors.qualitative.Dark24,
                           ).data[1])
            project_allCells.write_html(self.project + '_' + step + '_' + 'allCells.html')
            # project_allCells.write_image(project + '_allCells.png', width=400, height=400)

            project_cumuCells = px.scatter(
                self.cumulative_results_m[self.cumulative_results_m['project'] == self.project], x='e', y=step,
                title=self.project + 'single point - receptor type average',
                color='type', template='presentation', width=800, height=800,
                marginal_x='rug', marginal_y='rug',
                color_discrete_sequence=px.colors.qualitative.Dark24,
            )
            project_cumuCells.add_trace(
                px.scatter(self.cumulative_results_m[self.cumulative_results_m['project'] == self.project], x='e',
                           y=step,
                           trendline='ols',
                           color_discrete_sequence=px.colors.qualitative.Dark24,
                           ).data[1])
            project_cumuCells.write_html(self.project + '_' + step + '_' + 'cumuCells.html')
            # project_cumuCells.write_image(self.project + '_cumuCells.png', width=400, height=400)

    def update_config(self):

        config = pd.read_csv(self.config)
        new_starters = self.cumulative_results_m.drop(columns=['project', 'model'])

        new_config = pd.merge(new_starters, config.loc[:, ['type', 'file', 'file_scn', 'tres', 'tcrit', 'model']],
                              on='type')
        new_config.to_csv('hjcfit_config_' + self.project + '_iter.csv')


parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", type=str)
args = parser.parse_args()

project = '_'.join(re.split('[_,.]', args.config)[2:-1])
ratesFile = 'hjcfit_rates_' + project + '.csv'

# project = 'bambi_test'
# ratesFile = 'hjcfit_rates_f64_bambi_test.csv'                                                                           # test file

ploter = REFERPloter(ratesFile, project, args.config)
ploter.REFER_plot_Auerbach()
ploter.update_config()
# ploter.REFER_plot()
