import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import statsmodels
import argparse


class REFERPloter:

    def __init__(self, rates_tabel, order):

        self.table_results = pd.read_csv(rates_tabel).drop(columns='Unnamed: 0')

        # HOWTO: calculate mean values in given categories
        # 1. pivot by values to aggregate, gives multiindex series with single column of values from 'values'
        # 2. index reset I: unstack by level=0, being the index level named by 'values'
        # 3. index reset II: reset remaining, unpivoted levels of index, those specified as 'columns' during pivot

        self.cumulative_results = pd.pivot_table(self.table_results, columns=['type', 'project'], values=['forward', 'equilibrium'], aggfunc=np.mean)
        self.cumulative_results = self.cumulative_results.unstack(level=0)
        self.cumulative_results.reset_index(inplace=True)

        self.table_results.sort_values('type', ascending=False, inplace=True)
        print(self.table_results)
        self.cumulative_results.sort_values('type', ascending=False, inplace=True)
        print(self.cumulative_results)

    def REFER_plot(self):

        # for each project separately

        for project in args.projectsOrder:

            project_allCells = px.scatter(self.table_results[self.table_results['project'] == project], x='equilibrium', y='forward',
                                          title='single point - single cell',
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
                                           title='single point - receptor type average',
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
                                          category_orders={'project': args.projectsOrder},
                                          color_discrete_sequence=px.colors.qualitative.Dark24,
                                          )
        allProjects_allCells.write_html('allProjects_allCells.html')
        allProjects_allCells.write_image('allProjects_allCells.png')

        allProjects_cumuCells = px.scatter(self.cumulative_results, x='equilibrium', y='forward',
                                          title='single point - receptor type average',
                                          color='project', trendline='ols', template='presentation', width=800, height=800,
                                          hover_data=['type'],
                                          category_orders={'project': args.projectsOrder},
                                          color_discrete_sequence=px.colors.qualitative.Dark24,
                                          )
        allProjects_cumuCells.write_html('allProjects_cumuCells.html')
        allProjects_cumuCells.write_image('allProjects_cumuCells.png')


# 'F200', 'F64', 'F45', 'F14', 'F31', 'H55', 'G254', 'F14F31'

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--ratesFile", type=str)
parser.add_argument("-c", "--configFile", type=str)
parser.add_argument("-p", "--projectsOrder", type=str, nargs='+')
args = parser.parse_args()

ploter = REFERPloter(args.ratesFile, args.projectsOrder, args.config)
ploter.REFER_plot()
