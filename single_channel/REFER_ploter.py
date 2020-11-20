import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import statsmodels
import argparse
import re


class REFERPloter:

    def __init__(self, rates_tabel, order, config):

        self.config = config
        self.order = order

        self.table_results = pd.read_csv(rates_tabel).drop(columns='Unnamed: 0')

        # WRONG !!!
        # HOWTO: calculate mean values in given categories WRONG
        # 1. pivot by values to aggregate, gives multiindex series with single column of values from 'values'
        # 2. index reset I: unstack by level=0, being the index level named by 'values'
        # 3. index reset II: reset remaining, unpivoted levels of index, those specified as 'columns' during pivot

        self.cumulative_results = self.table_results.drop(columns=['file'])
        # self.cumulative_results = pd.pivot_table(self.cumulative_results, columns=['project', 'type', 'model'], aggfunc=np.mean)
        # self.cumulative_results = self.cumulative_results.unstack(level=0)
        # self.cumulative_results = self.cumulative_results.reset_index(inplace=False)

        self.cumulative_results = self.cumulative_results.groupby(['project', 'type', 'model'], as_index=False)['alpha', 'beta', 'gamma', 'delta', 'equilibrium', 'forward'].mean()

        self.table_results.sort_values('type', ascending=False, inplace=True)
        self.cumulative_results.sort_values('type', ascending=False, inplace=True)

        print(self.table_results)
        print(self.cumulative_results)


    def REFER_plot(self):

        # for each project separately

        for project in self.order:

            project_allCells = px.scatter(self.table_results[self.table_results['project'] == project], x='equilibrium', y='forward',
                                          title=project + 'single point - single cell',
                                          color='type', template='presentation', width=800, height=800,
                                          marginal_x='rug', marginal_y='rug',
                                          color_discrete_sequence=px.colors.qualitative.Dark24,
                                          )
            project_allCells.add_trace(px.scatter(self.table_results[self.table_results['project'] == project], x='equilibrium', y='forward',
                                                  trendline='ols',
                                                  color_discrete_sequence=px.colors.qualitative.Dark24,
                                                  ).data[1])
            project_allCells.write_html(project + '_allCells.html')
            project_allCells.write_image(project + '_allCells.png')

            project_cumuCells = px.scatter(self.cumulative_results[self.cumulative_results['project'] == project], x='equilibrium', y='forward',
                                           title=project + 'single point - receptor type average',
                                           color='type', template='presentation', width=800, height=800,
                                           marginal_x='rug', marginal_y='rug',
                                           color_discrete_sequence=px.colors.qualitative.Dark24,
                                           )
            project_cumuCells.add_trace(px.scatter(self.cumulative_results[self.cumulative_results['project'] == project], x='equilibrium', y='forward',
                                                   trendline='ols',
                                                   color_discrete_sequence=px.colors.qualitative.Dark24,
                                                   ).data[1])
            project_cumuCells.write_html(project + '_cumuCells.html')
            project_cumuCells.write_image(project + '_cumuCells.png')

        # for all projects

        allProjects_allCells = px.scatter(self.table_results, x='equilibrium', y='forward',
                                          title='single point - single cell',
                                          color='project', trendline='ols', template='presentation', width=800, height=800,
                                          hover_data=['type'],
                                          category_orders={'project': self.order},
                                          color_discrete_sequence=px.colors.qualitative.Dark24,
                                          )
        allProjects_allCells.write_html('allProjects_allCells.html')
        allProjects_allCells.write_image('allProjects_allCells.png')

        allProjects_cumuCells = px.scatter(self.cumulative_results, x='equilibrium', y='forward',
                                          title='single point - receptor type average',
                                          color='project', trendline='ols', template='presentation', width=800, height=800,
                                          hover_data=['type'],
                                          category_orders={'project': self.order},
                                          color_discrete_sequence=px.colors.qualitative.Dark24,
                                          )
        allProjects_cumuCells.write_html('allProjects_cumuCells.html')
        allProjects_cumuCells.write_image('allProjects_cumuCells.png')

    def update_config(self):

        # works only for single project mode

        config = pd.read_csv(self.config)
        new_starters = self.cumulative_results.drop(columns=['project', 'model', 'equilibrium', 'forward'])

        new_config = pd.merge(new_starters, config.loc[:, ['type', 'file', 'file_scn', 'tres', 'tcrit', 'model']], on='type')
        new_config.to_csv('hjcfit_config_' + self.order[0] + '_iter.csv')

# 'F200', 'F64', 'F45', 'F14', 'F31', 'H55', 'G254', 'F14F31'


parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", type=str)
parser.add_argument("-p", "--projectsOrder", type=str, nargs='+')
args = parser.parse_args()

project = '_'.join(re.split('[_,.]', args.config)[2:-1])
ratesFile = 'hjcfit_rates_' + project + '.csv'
if args.projectsOrder is not None:
    projectsOrder = args.projectsOrder
else:
    projectsOrder = [project]

ploter = REFERPloter(ratesFile, projectsOrder, args.config)
ploter.update_config()
ploter.REFER_plot()
