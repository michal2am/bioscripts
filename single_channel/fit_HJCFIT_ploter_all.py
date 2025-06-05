import numpy as np
import pandas as pd
import plotly.express as px
import statsmodels
import argparse
import re
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import math
import matplotlib.ticker as ticker


parser = argparse.ArgumentParser()
parser.add_argument("-dfs", "--data_files", nargs='+')
parser.add_argument("--control")
args = parser.parse_args()

datas_hjcfit= pd.concat([pd.read_csv(df, index_col=0) for df in args.data_files])
# print(datas_hjcfit)
# below may be used to truncate project name to residue id only:
# datas_hjcfit['project'] = datas_hjcfit.apply(lambda cell: cell['project'].split('_')[0], axis=1)

# print(datas_hjcfit[datas_hjcfit['type'] == 'WT'][['project', 'type', 'model', 'alpha', 'beta']])


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
                           labels={'equilibrium_raw': 'log(equilibrium constant)',
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


def plot_each_project_single_wt_seaborn(wt):

    WT_project = wt

    forward = 'new_forward'
    equilibrium = 'new_equilibrium'

    phis = []
    projects = []
    wts = []

    for project in datas_hjcfit['project'].unique():
        single_project = datas_hjcfit[(datas_hjcfit['project'] == project) & ~(datas_hjcfit['type'] == 'WT')]
        selected_WT = datas_hjcfit[(datas_hjcfit['project'] == WT_project) & (datas_hjcfit['type'] == 'WT')]
        new_set = pd.concat([selected_WT, single_project])
        print(single_project)
        print(selected_WT)
        new_wt_a = selected_WT.alpha.values[0]
        new_wt_b = selected_WT.beta.values[0]
        new_set['new_forward'] = new_set.apply(lambda row: np.log(row['beta']/new_wt_b), axis=1)
        new_set['new_equilibrium'] = new_set.apply(lambda row: np.log((row['beta']/new_wt_b)/(row['alpha']/new_wt_a)), axis=1)

        slope, intercept, r_value, p_value, std_err= stats.linregress(new_set[equilibrium], new_set[forward])
        print("Slope: {}, intercept: {}, r_value: {}, p_value: {}, std_err: {}".format(slope, intercept, r_value, p_value, std_err))
        phis.append(slope)
        # print(single_project)
        projects.append(project)
        wts.append(wt)

        sns.set_style()
        sns.set_context("paper")

        g = sns.lmplot(
            data=new_set,
            x=equilibrium, y=forward,
            ci=None, #68,
            scatter=False,
            line_kws={"color": "black", 'zorder': 1},
            # palette=sns.xkcd_palette(["pale red", 'windows blue']),
            height=2, aspect=1,
        )

        g.map(sns.scatterplot, data=new_set, x=equilibrium, y=forward, hue='type', s=50)
        g.add_legend(bbox_to_anchor=(1, 0.5), loc='center right')

        g.set(#xlabel='log(equilibrium constant)',
              #ylabel='log(forward rate)',
              xlabel=None,
              ylabel=None,
              #xticks=np.linspace(math.floor(new_set["equilibrium"].min()), new_set["equilibrium"].max(), 3),
              #yticks=np.linspace(math.floor(new_set["forward"].min()), new_set["forward"].max(), 3))
              yticks=[-0.6, -0.3, 0.0],
              xticks=[-0.75, -0.375, 0.0])
        g.fig.suptitle("{:.2f} +/- {:.2f}".format(slope, std_err))
        #g.fig.subplots_adjust(top=0.8)

        ax = g.axes[0, 0]
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

        g.despine(trim=False)
        #g.legend.set_title("")

        plt.tight_layout(rect=[0, 0, 0.8, 1])

        plt.savefig('project_' + project + '_wt_' + WT_project  + '_sns.png', dpi=150)
        #plt.show()

    all_phis = pd.DataFrame({'project': projects, 'phi': phis, 'wt': wts})
    print(all_phis)
    all_phis.to_csv('phis_wt_' + wt + '.csv')

def plot_each_project_no_wt_seaborn():

    forward = 'forward_raw'
    equilibrium = 'equilibrium_raw'
    phis = []
    projects = []
    wts = []

    for project in datas_hjcfit['project'].unique():
        single_project = datas_hjcfit[(datas_hjcfit['project'] == project) & ~(datas_hjcfit['type'] == 'WT')]

        slope, intercept, r_value, p_value, std_err = stats.linregress(single_project[equilibrium], single_project[forward])
        print(slope)
        phis.append(slope)
        print(single_project)
        projects.append(project)
        wts.append('wt_no')

        sns.set_style()
        sns.set_context("talk")

        g = sns.lmplot(
            data=single_project,
            x=equilibrium, y=forward,
            ci=68,
            scatter=False,
            # palette=sns.xkcd_palette(["pale red", 'windows blue']),
            height=4, aspect=1,
        )

        g.map(sns.scatterplot, data=single_project, x=equilibrium, y=forward, hue='type')
        g.add_legend()
        g.set(xlabel='log(equilibrium constant)',
              ylabel='log(forward rate)')
        g.fig.suptitle("{}: {:.2f}".format(project.upper(), slope))
        g.fig.subplots_adjust(top=0.8)
        g.despine(trim=False)
        #g.legend.set_title("")

        plt.savefig('project_' + project + '_wt_no' + '_sns.png', dpi=300)
        #plt.show()

    all_phis = pd.DataFrame({'project': projects, 'phi': phis, 'wt': wts})
    print(all_phis)
    all_phis.to_csv('phis_wt_no.csv')

plot_each_project_single_wt_seaborn(args.control)
# plot_each_project_single_wt_seaborn('f14')
# plot_each_project_no_wt_seaborn()