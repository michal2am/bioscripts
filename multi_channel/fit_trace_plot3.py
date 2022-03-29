import pyabf
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import plotly.express as px

import pandas as pd
import seaborn as sns


sim_30_traces = ['WT_fl.abf', 'P277K_fl.abf', 'P277H_fl.abf', 'P277H_poprawki_fl.abf']                                  # those should be abf 2.0+ files
sim_30_titles = ['P277_WT', 'P277K', 'P277H', 'P277H_pop']
sim_30_colors = [(0, 0, 0), (0.88, 0.54, 0.31), (0.8, 0.74, 0.43), 'black']

#sim_30_LC_traces = ['WT_LC_30_final.abf', 'G254V_LC_30_final.abf']
#sim_30_LC_titles = ['WT_LC', 'G254V']
#sim_30_LC_colors = [(0, 0, 0), (0.77, 0.20, 0.25)]

#sim_500_traces = ['WT_500_final.abf', 'L300V_500_final.abf']
#sim_500_titles = ['WT_500', 'L300V']
#sim_500_colors = [(0, 0, 0), (0.60, 0.75, 0.75)]

sns.set_style()
sns.set_context("paper")

def plot_traces(sim_traces, sim_titles, sim_colors, window):

    traces = pd.DataFrame()

    for trace, type in zip(sim_traces, sim_titles):

        print('Reading file {}'.format(trace))
        sim_file = pyabf.ABF(trace)


        #exp_file.setSweep(sweepNumber)

        sim_t = sim_file.sweepX
        sim_p = sim_file.sweepY / max(sim_file.sweepY)

        trace_pd = pd.DataFrame()
        trace_pd['t'] = sim_t
        trace_pd['p'] = sim_p
        trace_pd['type'] = type

        traces = traces.append(trace_pd)

    trace_30_window = traces[(traces['t'] < window[1]) & (traces['t'] > window[0])]
    trace_30_window.to_csv('test.csv')

    fig = px.line(trace_30_window, x='t', y='p', color='type')
    fig.show()

    sns.relplot(x='t', y='p', hue='type', kind="line", data=trace_30_window,
                palette=sim_colors,
                height=1.54, aspect=2.,
                )
    #sns.relplot(x='t', y='p', kind="line", data=exp_traces)
    plt.savefig(sim_titles[0] + '.png', dpi=600)
    plt.show()

plot_traces(sim_30_traces, sim_30_titles, sim_30_colors, (0.04, 0.08))
#plot_traces(sim_30_LC_traces, sim_30_LC_titles, sim_30_LC_colors, (0.04, 0.08))
#plot_traces(sim_500_traces, sim_500_titles, sim_500_colors, (0, 1))