import pyabf
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import pandas as pd
import seaborn as sns

traces = pd.DataFrame()

sim = pyabf.ATF('fit_trace.atf')
sim_t = sim.sweepX
sim_p = sim.sweepY/max(sim.sweepY)

exp = pyabf.ATF('exp_trace.atf')
exp_t = exp.sweepX

exp_p = exp.sweepY

exp_start = [i for i, x in enumerate(exp_p > 0.01) if x][0]
sim_start = 497

delta_start = exp_start - sim_start

exp_p = list(exp_p[delta_start:]) + [0]*delta_start
exp_p /= max(exp_p)

traces['t'] = sim_t
traces['sim'] = sim_p
traces['exp'] = exp_p

traces = pd.melt(traces, id_vars=['t'])


print(traces)

sns.relplot(x='t', y='value', hue='variable', kind="line", data=traces);
plt.show()
