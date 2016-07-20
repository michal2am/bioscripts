import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import time as tm
from mm_solver_ode import SolverOde
from mm_solver_glp import SolverGlp
from mm_kinetic import Kinetic
from mm_kinetic import Stimulus

# model definition #

all_time = tm.time()

stimulus = Stimulus.pair_square([0., 7.], [5., 12.], 10.0)
suspend = [[5., 7.], [12., 20.]]

print("### New rate added:")

r_kon   = Kinetic.State.Rate('kon',     2.00, stimulus)
r_2kon  = Kinetic.State.Rate('2kon',    4.00, stimulus)
r_koff  = Kinetic.State.Rate('koff',    1.00)
r_2koff = Kinetic.State.Rate('2koff',   2.00)
r_d     = Kinetic.State.Rate('d',       1.00)
r_r     = Kinetic.State.Rate('r',       0.20)
r_b     = Kinetic.State.Rate('b',       3.00)
r_a     = Kinetic.State.Rate('a',       2.00)
r_0     = Kinetic.State.Rate('block',   0.00)

print("### New state added:")

st_r    = Kinetic.State(0, 'R',   [r_0,     r_2kon,   r_0,    r_0,  r_0],   1, False)
st_ar   = Kinetic.State(1, 'AR',  [r_koff,  r_0,      r_kon,  r_0,  r_0],   0, False)
st_a2r  = Kinetic.State(2, 'A2R', [r_0,     r_2koff,  r_0,    r_d,  r_b],   0, False)
st_a2d  = Kinetic.State(3, 'A2D', [r_0,     r_0,      r_r,    r_0,  r_0],   0, False)
st_a2o  = Kinetic.State(4, 'A2O', [r_0,     r_0,      r_a,    r_0,  r_0],   0, True)

jwm = Kinetic([st_r, st_ar, st_a2r, st_a2d, st_a2o])
states_na = [state.name for state in [st_r, st_ar, st_a2r, st_a2d, st_a2o]]


# solving!

solve_ode = True
solve_glp = True
solve_sto = False
plot_stim = True

ini_conc = np.array([[1.0, 0.0, 0.0, 0.0, 0.0]])
states_no = ini_conc.shape[1]
part = 100
t0 = 0
te = 20

# figure preparation

sns.set_style("ticks", {'legend.frameon': True})
sns.set_context("talk")
sns.set_palette('deep', states_no)
colors = sns.color_palette('deep', states_no)

fig1, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, gridspec_kw={'height_ratios': [7, 1]}, )
fig1.set_size_inches(10, 8)

if solve_ode:

    solver = SolverOde(jwm.trmn, ini_conc, t0, te)
    T, P = solver.get_results()

    T = np.array([T]).transpose()
    conc = np.concatenate((T, P), axis=1)
    pandas = pd.DataFrame(data=conc, columns=['time'] + states_na)
    pandas.plot(x='time', ax=ax1)


if solve_glp:
    solver = SolverGlp(jwm, ini_conc, part, t0, te, suspend)

    allT, allP, mcT, mcP = solver.get_results()

    P = allP[0]
    T = allT[0]
    Pn = np.array([[float('nan') if state == 0.0 else state for state in step] for step in P])
    T = np.array([T]).transpose()

    # all state dot plot
    conc = np.concatenate((T, Pn), axis=1)
    pandas = pd.DataFrame(data=conc, columns=['time'] + states_na)
    pandas.plot(x='time', style=['o']*5, ax=ax1)

    # selected state step plot
    Pnplot_no = 4
    Pnplot_na = 'A2O'
    conc = np.concatenate((T, P), axis=1)
    pandas = pd.DataFrame(data=conc, columns=['time'] + states_na)
    pandas.plot(x='time', y=Pnplot_na, drawstyle='steps-post', color=colors[Pnplot_no], ax=ax1, linewidth=1.5)

    # cumulative trajectories
    mcT = np.array([mcT]).transpose()
    conc = np.concatenate((mcT, mcP), axis=1)
    pandas = pd.DataFrame(data=conc, columns=['time'] + states_na)
    pandas.plot(x='time', ax=ax1)

if solve_sto:
    pass

if plot_stim:
    T = np.linspace(t0 - 0.25, te, 1e5)
    stim = [stimulus(ti) for ti in T]
    ax2.step(T, stim, 'k-')

    ax2.set_xlim([t0 - 0.5, te + 0.5])
    ax2.set_ylim([0, max(stim) + 0.1])
    ax2.get_yaxis().set_visible(False)
    ax2.set_xlabel('time [ms]')
    ax2.set_ylabel('stimulus')

ax1.set_xlim([t0 - 0.5, te + 0.5])
ax1.set_ylim([0, 1.05])
ax1.legend(states_na)
ax1.set_ylabel('state probability')
sns.despine()
plt.show()

print("--- %s seconds ---" % (tm.time() - all_time))
