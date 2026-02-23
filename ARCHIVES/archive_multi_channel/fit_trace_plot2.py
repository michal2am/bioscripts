import pyabf
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import pandas as pd
import seaborn as sns

traces = pd.DataFrame()

'''
for sim in ['11', '12b', '24', '8', '12a', '13', '6', '9']:

    sim_file = pyabf.ABF(sim + '_sim.abf')
    sim_t = sim.sweepX
    sim_p = sim.sweepY/max(sim.sweepY)

    trace_pd = pd.DataFrame()
    trace_pd['t'] = sim_t
    trace_pd['p'] = sim_p
    trace_pd['type'] = 'sim'
    trace_pd['cell'] = sim

    print(trace_pd)
    
'''

exp_file = pyabf.ABF('all_exp.abf')
exp_traces = pd.DataFrame()

for sweepNumber in exp_file.sweepList:
    exp_file.setSweep(sweepNumber)

    exp_t = exp_file.sweepX
    exp_p = exp_file.sweepY/max(exp_file.sweepY)

    trace_pd = pd.DataFrame()
    trace_pd['t'] = exp_t
    trace_pd['p'] = exp_p
    trace_pd['type'] = 'exp'
    trace_pd['cell'] = sweepNumber

    print(trace_pd)
    exp_traces = exp_traces.append(trace_pd)

print(exp_traces)

sns.relplot(x='t', y='p', hue='cell', kind="line", data=exp_traces)
#sns.relplot(x='t', y='p', kind="line", data=exp_traces)
plt.show()
