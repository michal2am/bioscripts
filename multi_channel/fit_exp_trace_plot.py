import pyabf
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import pandas as pd
import seaborn as sns


traces = ['8_sim_1d_n.abf', '2018_05_16_0008_reduced_1_toModel_n.abf']
titles = ['sim', 'exp']
colors = ['black', 'red']


sns.set_style()
sns.set_context("paper")

def plot_traces(sim_traces, sim_titles, sim_colors, window):

    traces = pd.DataFrame()

    for trace, type in zip(sim_traces, sim_titles):

        sim_file = pyabf.ABF(trace)


        #exp_file.setSweep(sweepNumber)

        sim_t = sim_file.sweepX
        sim_p = sim_file.sweepY / max(sim_file.sweepY)

        trace_pd = pd.DataFrame()
        trace_pd['t'] = sim_t
        trace_pd['p'] = sim_p
        trace_pd['type'] = type

        traces = traces.append(trace_pd)

    print(traces)
    trace_30_window = traces[(traces['t'] < window[1]) & (traces['t'] > window[0])]

    sns.relplot(x='t', y='p', hue='type', kind="line", data=trace_30_window,
                palette=sim_colors,
                height=1.54, aspect=2.,
                )
    #sns.relplot(x='t', y='p', kind="line", data=exp_traces)
    plt.savefig(sim_titles[0] + '.png', dpi=600)
    plt.show()

plot_traces(traces, titles, colors, (0.0, 2.0))