import numpy as np
import pandas as pd
import plotly.express as px
import statsmodels
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-dfs", "--data_files", nargs='+')
args = parser.parse_args()

datas_hjcfit= pd.concat([pd.read_csv(df, index_col=0) for df in args.data_files])
datas_hjcfit['project'] = datas_hjcfit.apply(lambda cell: cell['project'].split('_')[0], axis=1)

print(datas_hjcfit[datas_hjcfit['type'] == 'WT'])


def plot_each_project_single_wt(wt):

    WT_project = wt

    for project in datas_hjcfit['project'].unique():
        single_project = datas_hjcfit[(datas_hjcfit['project'] == project) & ~(datas_hjcfit['type'] == 'WT')]
        selected_WT = datas_hjcfit[(datas_hjcfit['project'] == WT_project) & (datas_hjcfit['type'] == 'WT')]
        new_set = pd.concat([selected_WT, single_project])
        print(new_set)

        plot = (px.scatter(new_set, x='equilibrium_raw', y='forward_raw',
                           # title="{} REFER by HJCFIT CO, WT from {}".format(project, WT_project),
                           title="{}".format(project.upper()),
                           labels={'equilibrium_raw': 'log(equilibrium rate)',
                                   'forward_raw': 'log(forward rate)'},
                           color='type', template='presentation', width=400, height=400,
                           hover_name='type', hover_data=['alpha', 'beta'],
                           color_discrete_sequence=px.colors.qualitative.Dark24,
                           ))
        plot.add_trace(
            px.scatter(new_set, x='equilibrium_raw', y='forward_raw',
                       trendline='ols',
                       color_discrete_sequence=px.colors.qualitative.Dark24,
                       ).data[1])
        #plot.show()

        plot.write_html('project_' + project + '_wt_' + WT_project  + '_co.html')
        plot.write_image('project_' + project + '_wt_' + WT_project  + '.png', width=400, height=400, scale=1)

plot_each_project_single_wt('e153')