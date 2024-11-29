# plotowanie przebiegu eksperymentalnego z fitem channellaba
# pliki wejściowe to sklejone kolumny plików .fit

import pyabf
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import pandas as pd
import seaborn as sns


#sim_fit_5_traces = pd.read_csv('traces_5.csv', sep=',')
#sim_fit_5_traces = sim_fit_5_traces.melt(id_vars='t', var_name='type', value_name='popen')
#print(sim_fit_5_traces)

sim_fit_30_traces = pd.read_csv('fitted10.csv', sep=',')
sim_fit_30_traces = sim_fit_30_traces.melt(id_vars='t', var_name='type', value_name='popen')
print(sim_fit_30_traces)

#sim_fit_1000_traces = pd.read_csv('traces_1000.csv')
#sim_fit_1000_traces = sim_fit_1000_traces.melt(id_vars='t', var_name='type', value_name='popen')
#print(sim_fit_1000_traces)

'''
#VERY OLD#

sim_30_traces = ['WT_30_final.abf', 'G258V_30_final.abf', 'L296V_30_final.abf']
sim_30_titles = ['WT', 'G258V', 'L296V']
sim_30_colors = [(0, 0, 0), (0.88, 0.54, 0.31), (0.8, 0.74, 0.43)]

sim_30_LC_traces = ['WT_LC_30_final.abf', 'G254V_LC_30_final.abf']
sim_30_LC_titles = ['WT_LC', 'G254V']
sim_30_LC_colors = [(0, 0, 0), (0.77, 0.20, 0.25)]

sim_500_traces = ['WT_500_final.abf', 'L300V_500_final.abf']
sim_500_titles = ['WT_500', 'L300V']
sim_500_colors = [(0, 0, 0), (0.60, 0.75, 0.75)]
'''

sns.set_style()
sns.set_context("paper")

def plot_traces5(sim_fit_5_traces, sim_fit_colors, simt_fit_types, file):

    for_plot = sim_fit_5_traces.loc[sim_fit_5_traces['type'].isin(simt_fit_types)]

    sns.relplot(x='t', y='popen', hue='type', kind="line", data=for_plot,
                palette=sim_fit_colors,
                #style='type',
                height=2, aspect=1.,
                #linewidth=3,
                #scale=2,
                size='type',
                sizes=(0.5, 1)
                )
    plt.savefig(file + '.png', dpi=600)
    #plt.show()

def plot_traces30(sim_fit_30_traces, sim_fit_colors, simt_fit_types, file):

    for_plot = sim_fit_30_traces.loc[sim_fit_30_traces['type'].isin(simt_fit_types)]

    sns.relplot(x='t', y='popen', hue='type', kind="line", data=for_plot,
                palette=sim_fit_colors,
                #style='type',
                height=2, aspect=2,
                #linewidth=6,
                #scale=2,
                size='type',
                sizes=(0.5, 2)
                )
    plt.savefig(file + '.png', dpi=600)
    #plt.show()

def plot_traces1000(sim_fit_1000_traces, sim_fit_colors, simt_fit_types, file):

    for_plot = sim_fit_1000_traces.loc[sim_fit_1000_traces['type'].isin(simt_fit_types)]

    sns.relplot(x='t', y='popen', hue='type', kind="line", data=for_plot,
                palette=sim_fit_colors,
                #style='type',
                height=2, aspect=2.,
                #linewidth=3,
                size='type',
                sizes=(0.5, 1)
                )
    plt.savefig(file + '.png', dpi=600)
    #plt.show()


#plot_traces5(sim_fit_5_traces, [(0.57, 0.55, 0.55), 'black'], ['WT_2D_exp', 'WT_2D_fit'], 'WT_2D')
#plot_traces5(sim_fit_5_traces, [(0.40, 0.67, 0.67), 'black'], ['L300V_exp', 'L300V_fit'], 'L300V')


plot_traces30(sim_fit_30_traces, [(0.57, 0.55, 0.55), 'black'], ['WT_EXP', 'WT_FIT'], 'WT')
#plot_traces30(sim_fit_30_traces, [(0.57, 0.55, 0.55), 'black'], ['WT_exp2', 'WT_fit2'], 'WT2')
#plot_traces30(sim_fit_30_traces, [(0.57, 0.55, 0.55), 'black'], ['WT_LC_exp', 'WT_LC_fit'], 'WT_LC')
#plot_traces30(sim_fit_30_traces, [(0.88, 0.54, 0.31), 'black'], ['k_exp', 'k_fit'], 'P277K')
#plot_traces30(sim_fit_30_traces, [(0.70, 0.62, 0.17), 'black'], ['h_exp', 'h_fit'], 'P277H')
#plot_traces30(sim_fit_30_traces, [(0.77, 0.20, 0.25), 'black'], ['G254V_exp', 'G254V_fit'], 'G254V')


#plot_traces1000(sim_fit_1000_traces, ['grey', 'black'], ['WT_exp', 'WT_fit'], 'WT_1000')
#plot_traces1000(sim_fit_1000_traces, [(0.40, 0.67, 0.67), 'black'], ['L300V_exp', 'L300V_fit'], 'L300V_1000')


#plot_traces(sim_30_LC_traces, sim_30_LC_titles, sim_30_LC_colors, (0.04, 0.08))
#plot_traces(sim_500_traces, sim_500_titles, sim_500_colors, (0, 1))