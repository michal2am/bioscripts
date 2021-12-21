import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import statsmodels
import argparse
import re
from plotly.subplots import make_subplots



class REFERPloter:

    def __init__(self, rates_table, project, config):

        self.config = config
        self.project = project

        self.table_results = pd.read_csv(rates_table).drop(columns='Unnamed: 0')
        self.wt_results = pd.read_csv('hjcfit_rates_WT_pooled.csv').drop(columns='Unnamed: 0')

        print(self.table_results)
        print(self.wt_results)

        self.table_results = pd.concat([self.table_results, self.wt_results])
        self.table_results['project'] = self.project

        print(self.table_results)

        self.cumulative_results = self.table_results.drop(columns=['file'])

        # means
        self.cumulative_results_m = self.cumulative_results.groupby(['project', 'type', 'model'], as_index=False)['alpha', 'beta', 'gamma', 'delta', 'e', 'f1', 'f2', 'f3'].mean()
        # % errors
        self.cumulative_results_d = self.cumulative_results.groupby(['project', 'type', 'model'], as_index=False)['alpha', 'beta', 'gamma', 'delta', 'e', 'f1', 'f2', 'f3'].apply(lambda x: np.std(x)/np.mean(x)*100)

        self.table_results.sort_values('type', ascending=False, inplace=True)
        self.cumulative_results_m.sort_values('type', ascending=False, inplace=True)
        self.cumulative_results_d.sort_values('type', ascending=False, inplace=True)

        print(self.table_results)
        print(self.cumulative_results_m)
        print(self.cumulative_results_d)

        self.cumulative_results_m.to_csv('hjcfit_rates_m_' + self.project + '.csv')
        self.cumulative_results_d.to_csv('hjcfit_rates_d_' + self.project + '.csv')

    def REFER_plot(self):


        for step in ['f1', 'f2', 'f3']:

                project_allCells = (px.scatter(self.table_results[self.table_results['project'] == self.project], x='e', y=step,
                                              title=self.project + 'single point - single cell',
                                              color='type', template='presentation', width=800, height=800,
                                              marginal_x='rug', marginal_y='rug',
                                              hover_name='file', hover_data=['alpha', 'beta', 'gamma', 'delta'],
                                              color_discrete_sequence=px.colors.qualitative.Dark24,
                                              ))
                project_allCells.add_trace(px.scatter(self.table_results[self.table_results['project'] == self.project], x='e', y=step,
                                                      trendline='ols',
                                                      color_discrete_sequence=px.colors.qualitative.Dark24,
                                                      ).data[1])
                project_allCells.write_html(self.project + '_' + step + '_' + 'allCells.html')
                # project_allCells.write_image(project + '_allCells.png', width=400, height=400)

                project_cumuCells = px.scatter(self.cumulative_results_m[self.cumulative_results_m['project'] == self.project], x='e', y=step,
                                               title=self.project + 'single point - receptor type average',
                                               color='type', template='presentation', width=800, height=800,
                                               marginal_x='rug', marginal_y='rug',
                                               color_discrete_sequence=px.colors.qualitative.Dark24,
                                               )
                project_cumuCells.add_trace(px.scatter(self.cumulative_results_m[self.cumulative_results_m['project'] == self.project], x='e', y=step,
                                                       trendline='ols',
                                                       color_discrete_sequence=px.colors.qualitative.Dark24,
                                                       ).data[1])
                project_cumuCells.write_html(self.project + '_' + step + '_' + 'cumuCells.html')
                # project_cumuCells.write_image(self.project + '_cumuCells.png', width=400, height=400)

    def update_config(self):

        config = pd.read_csv(self.config)
        new_starters = self.cumulative_results_m.drop(columns=['project', 'model', 'e', 'f1', 'f2', 'f3'])

        new_config = pd.merge(new_starters, config.loc[:, ['type', 'file', 'file_scn', 'tres', 'tcrit', 'model']], on='type')
        new_config.to_csv('hjcfit_config_' + self.project + '_iter.csv')


parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", type=str)
args = parser.parse_args()

project = '_'.join(re.split('[_,.]', args.config)[2:-1])
ratesFile = 'hjcfit_rates_' + project + '.csv'


ploter = REFERPloter(ratesFile, project, args.config)
ploter.update_config()
ploter.REFER_plot()
