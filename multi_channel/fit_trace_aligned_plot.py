import pyabf
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import pandas as pd
import seaborn as sns


sim_fit_30_traces = pd.read_csv('traces_30.csv')
sim_fit_30_traces = sim_fit_30_traces.melt(id_vars='t', var_name='type', value_name='popen')
print(sim_fit_30_traces)

sim_fit_1000_traces = pd.read_csv('traces_1000.csv')
sim_fit_1000_traces = sim_fit_1000_traces.melt(id_vars='t', var_name='type', value_name='popen')
print(sim_fit_1000_traces)

'''
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

def plot_traces30(sim_fit_30_traces, sim_fit_colors, simt_fit_types, file):

    for_plot = sim_fit_30_traces.loc[sim_fit_30_traces['type'].isin(simt_fit_types)]

    sns.relplot(x='t', y='popen', hue='type', kind="line", data=for_plot,
                palette=sim_fit_colors,
                style='type',
                height=2, aspect=1.,
                size='type',
                )
    plt.savefig(file + '.png', dpi=600)
    plt.show()

def plot_traces1000(sim_fit_30_traces, sim_fit_colors, simt_fit_types, file):

    for_plot = sim_fit_30_traces.loc[sim_fit_30_traces['type'].isin(simt_fit_types)]

    sns.relplot(x='t', y='popen', hue='type', kind="line", data=for_plot,
                palette=sim_fit_colors,
                style='type',
                height=2, aspect=2.,
                size='type',
                )
    plt.savefig(file + '.png', dpi=600)
    plt.show()


plot_traces30(sim_fit_30_traces, ['grey', 'black'], ['WT_exp', 'WT_fit'], 'WT')
plot_traces30(sim_fit_30_traces, ['grey', 'black'], ['WT_LC_exp', 'WT_LC_fit'], 'WT_LC')

plot_traces30(sim_fit_30_traces, [(0.88, 0.54, 0.31), 'black'], ['G258V_exp', 'G258V_fit'], 'G258V')
plot_traces30(sim_fit_30_traces, [(0.8, 0.74, 0.43), 'black'], ['L296V_exp', 'L296V_fit'], 'L296V')
plot_traces30(sim_fit_30_traces, [(0.77, 0.20, 0.25), 'black'], ['G254V_exp', 'G254V_fit'], 'G254V')
plot_traces30(sim_fit_30_traces, [(0.60, 0.75, 0.75), 'black'], ['L300V_exp', 'L300V_fit'], 'L300V')


plot_traces1000(sim_fit_1000_traces, ['grey', 'black'], ['WT_exp', 'WT_fit'], 'WT_1000')
plot_traces1000(sim_fit_1000_traces, [(0.60, 0.75, 0.75), 'black'], ['L300V_exp', 'L300V_fit'], 'L300V_1000')


#plot_traces(sim_30_LC_traces, sim_30_LC_titles, sim_30_LC_colors, (0.04, 0.08))
#plot_traces(sim_500_traces, sim_500_titles, sim_500_colors, (0, 1))