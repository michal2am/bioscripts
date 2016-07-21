import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import seaborn as sns
import pandas as pd
import scipy.stats as sts
import time as tm
import logging as log
from mm_solver_ode import SolverOde
from mm_solver_glp import SolverGlp
from mm_kinetic import Kinetic
from mm_kinetic import Stimulus

log.basicConfig(filename='mm_kin.log', filemode='w', format='%(message)s', level=log.DEBUG)
log.info("### Kinetic solver starts!")

# model definition #

model = 'kisiel'
# model = 'jwm'

all_time = tm.time()

'''
stimulus = Stimulus.pair_square([0., 10.], [5., 15.], 10.0)
suspend = [[5., 10.], [15., 20.]]
'''

stimulus = Stimulus.square(0., 1.0, 0.1)
suspend = [[2.5], [10.]]


log.info("### Adding rates:")

if model == 'kisiel':

    r_kon   = Kinetic.State.Rate('kon',     2*1e7,   stimulus)
    r_2kon  = Kinetic.State.Rate('2kon',    2*2*1e7, stimulus)
    r_koff  = Kinetic.State.Rate('koff',    1116)
    r_2koff = Kinetic.State.Rate('2koff',   2*1116)
    r_b1    = Kinetic.State.Rate('b1',      15)
    r_b2    = Kinetic.State.Rate('b2',      200)
    r_b3    = Kinetic.State.Rate('b3',      35)
    r_b4    = Kinetic.State.Rate('b4',      200)
    r_b5    = Kinetic.State.Rate('b5',      100)
    r_a1    = Kinetic.State.Rate('a1',      20000)
    r_a2    = Kinetic.State.Rate('a2',      800)
    r_a3    = Kinetic.State.Rate('a3',      3333)
    r_a4    = Kinetic.State.Rate('a4',      17500)
    r_a5    = Kinetic.State.Rate('a5',      3000)
    r_d1    = Kinetic.State.Rate('d1',      310)
    r_d2    = Kinetic.State.Rate('d2',      800)
    r_d4    = Kinetic.State.Rate('d4',      1000)
    r_r1    = Kinetic.State.Rate('r1',      5)
    r_r2    = Kinetic.State.Rate('r2',      400)
    r_r4    = Kinetic.State.Rate('r4',      10)
    r_y2    = Kinetic.State.Rate('y2',      4400)
    r_g2    = Kinetic.State.Rate('g2',      4300)
    r_0     = Kinetic.State.Rate('block',   0.00)

if model == 'jwm':

    r_kon   = Kinetic.State.Rate('kon',     2.00, stimulus)
    r_2kon  = Kinetic.State.Rate('2kon',    4.00, stimulus)
    r_koff  = Kinetic.State.Rate('koff',    4.00)
    r_2koff = Kinetic.State.Rate('2koff',   8.00)
    r_d     = Kinetic.State.Rate('d',       1.00)
    r_r     = Kinetic.State.Rate('r',       0.70)
    r_b     = Kinetic.State.Rate('b',       5.00)
    r_a     = Kinetic.State.Rate('a',       3.00)
    r_0     = Kinetic.State.Rate('block',   0.00)

log.info("### Adding states:")

if model == 'kisiel':

    #                                  A1O      A2O     A3O     A40     A50     R       A1R     A2R     A2F     A1D     A2D     A4D
    st_a1o  = Kinetic.State(0, 'A1O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_a1,   r_0,    r_0,    r_d1,   r_0,    r_0],   0, True)
    st_a2o  = Kinetic.State(1, 'A2O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_a2,   r_0,    r_0,    r_0],   0, True)
    st_a3o  = Kinetic.State(2, 'A3O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_a3,   r_0,    r_0],   0, True)
    st_a4o  = Kinetic.State(3, 'A4O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_a4,   r_0,    r_0,    r_0,    r_d4],  0, True)
    st_a5o  = Kinetic.State(4, 'A5O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_a5],  0, True)
    st_r    = Kinetic.State(5, 'R ' , [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_2kon, r_0,    r_0,    r_0,    r_0,    r_0],   1, False)
    st_a1r  = Kinetic.State(6, 'A1R', [r_a1,    r_0,    r_0,    r_0,    r_0,    r_koff, r_0,    r_kon,  r_0,    r_0,    r_0,    r_0],   0, False)
    st_a2r  = Kinetic.State(7, 'A2R', [r_0,     r_0,    r_0,    r_b4,   r_0,    r_0,    r_2koff,r_0,    r_g2,   r_0,    r_0,    r_0],   0, False)
    st_a2f  = Kinetic.State(8, 'A2F', [r_0,     r_b2,   r_0,    r_0,    r_0,    r_0,    r_0,    r_y2,    r_0,   r_0,    r_d2,   r_0],   0, False)
    st_a1d  = Kinetic.State(9, 'A1D', [r_r1,    r_0,    r_b3,   r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0],   0, False)
    st_a2d  = Kinetic.State(10,'A2D', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_r2,   r_0,    r_0,    r_0],   0, False)
    st_a4d  = Kinetic.State(11,'A4D', [r_0,     r_0,    r_0,    r_r4,   r_b5,   r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0],   0, False)

if model == 'jwm':

    st_r    = Kinetic.State(0, 'R',   [r_0,     r_2kon,   r_0,    r_0,  r_0],   1, False)
    st_ar   = Kinetic.State(1, 'AR',  [r_koff,  r_0,      r_kon,  r_0,  r_0],   0, False)
    st_a2r  = Kinetic.State(2, 'A2R', [r_0,     r_2koff,  r_0,    r_d,  r_b],   0, False)
    st_a2d  = Kinetic.State(3, 'A2D', [r_0,     r_0,      r_r,    r_0,  r_0],   0, False)
    st_a2o  = Kinetic.State(4, 'A2O', [r_0,     r_0,      r_a,    r_0,  r_0],   0, True)


if model == 'kisiel':

    states = [st_a1o, st_a2o, st_a3o, st_a4o, st_a5o, st_r, st_a1r, st_a2r, st_a2f, st_a1d, st_a2d, st_a4d]
    model_kinetic = Kinetic(states)

if model == 'jwm':
    states = [st_r, st_ar, st_a2r, st_a2d, st_a2o]
    model_kinetic = Kinetic(states)


states_na = model_kinetic.states_na
states_no = model_kinetic.states_no
ini_conc = np.array([model_kinetic.ini_conc])
opsh = [model_kinetic.states_op, model_kinetic.states_sh]


# solving!

solve_ode = True
solve_glp = False
solve_sto = False
plot_stim = True


part = 10
t0 = 0
te = 2

# figure preparation

sns.set_style("ticks", {'legend.frameon': True})
sns.set_context("poster")
sns.set_palette('Paired', states_no)
colors = sns.color_palette('deep', states_no)

if solve_ode:

    fig_ode = plt.figure()
    fig_ode.canvas.set_window_title('ODE results')
    gs = grd.GridSpec(6, 2, width_ratios=[3, 2], hspace=2, wspace=0.3)
    ax2 = plt.subplot(gs[5, 0])
    ax1 = plt.subplot(gs[0:5, 0], sharex=ax2)
    ax4 = plt.subplot(gs[3:, 1])
    ax3 = plt.subplot(gs[0:3, 1])
    fig_ode.add_subplot(ax1, ax2, ax3, ax4)

    # ax1 trajectory
    solver = SolverOde(model_kinetic.trmn, ini_conc, states_na, t0, te, False)
    TP = solver.get_results()
    TP.plot(x='time', ax=ax1, legend=True)

    ax1.set_xlim([t0 - 0.5, te + 0.5])
    ax1.set_ylim([0, 1.00])
    ax1.set_ylabel('state probability')

    # ax2 stimulus
    T = np.linspace(t0 - 0.25, te, 1e5)
    stim = [stimulus(ti) for ti in T]
    ax2.step(T, stim, 'k-')

    ax2.set_xlim([t0 - 0.5, te + 0.5])
    ax2.set_ylim([0, max(stim)])
    ax2.set_yticks([0, max(stim)])

    ax2.set_ylabel('stimulus')
    ax2.set_xlabel('time [ms]')

    # ax3 distributions
    lifetimes = pd.DataFrame()
    time = np.linspace(0., .00001, num=1000)
    lifetimes['time'] = time
    sumrates = [np.sum(state_rates) for state_rates in model_kinetic.trm(0)]
    for no, name in enumerate(states_na):
        lifetimes[name] = sts.expon.pdf(time, loc=0, scale=1/sumrates[no])
    lifetimes.plot(x='time', ax=ax3, legend=False)

    ax3.set_ylabel('probability density function')
    #ax3.set_xlim([0, .31])
    ax3.set_xlabel('occupancy time [ms]')

    # ax4 equilibrium
    solver = SolverOde(model_kinetic.trmn, ini_conc, states_na, t0, te, True)
    TP = solver.get_results()[-1:].drop('time', 1)
    sns.barplot(data=TP, ax=ax4)

    ax4.set_ylabel('equilibrium occupancies')
    #ax4.set_xlim([0, 2.5])
    ax4.set_xlabel('state')

if solve_glp:
    solver = SolverGlp(model_kinetic, ini_conc, part, t0, te, suspend, opsh)

    allT, allP, mcT, mcP, distT, opshDistT = solver.get_results()

    P = allP[0]
    T = allT[0]
    Pn = np.array([[float('nan') if state == 0.0 else state for state in step] for step in P])
    T = np.array([T]).transpose()

    # all state dot plot
    conc = np.concatenate((T, Pn), axis=1)
    pandas = pd.DataFrame(data=conc, columns=['time'] + states_na)
    pandas.plot(x='time', style=['o']*5, ax=ax1, legend=False)

    # selected state step plot
    Pnplot_no = 4
    Pnplot_na = 'A2O'
    conc = np.concatenate((T, P), axis=1)
    pandas = pd.DataFrame(data=conc, columns=['time'] + states_na)
    pandas.plot(x='time', y=Pnplot_na, drawstyle='steps-post', color=colors[Pnplot_no], ax=ax1, linewidth=1.5, legend=False)

    # cumulative trajectories
    mcT = np.array([mcT]).transpose()
    conc = np.concatenate((mcT, mcP), axis=1)
    pandas = pd.DataFrame(data=conc, columns=['time'] + states_na)
    pandas.plot(x='time', ax=ax1, legend=False)

    # distributions
    for state_no, state in enumerate(states_na):
        loc, scale = sts.expon.fit(distT[state_no], floc=0)
        print(state, scale)
        x = np.linspace(0., 5., num=1000)
        y = sts.expon.pdf(x, loc=loc, scale=scale)
        pd.DataFrame({'x': x, state: y}).plot(x='x', ax=ax3)
        sns.distplot(distT[state_no], hist=False, rug=True, kde=False, ax=ax3, rug_kws={'height': 0.05})

    for state in opshDistT:
        loc, scale = sts.expon.fit(opshDistT[state], floc=0)
        print(state, scale)
        x = np.linspace(0., 5., num=1000)
        y = sts.expon.pdf(x, loc=loc, scale=scale)
        pd.DataFrame({'x': x, state: y}).plot(x='x', ax=ax4)
        sns.distplot(opshDistT[state], hist=False, rug=True, kde=False, ax=ax4, rug_kws={'height': 0.05})

if solve_sto:
    pass


#ax4.legend()

#ax3.legend(states_na)
sns.despine()

plt.show()

print("--- %s seconds ---" % (tm.time() - all_time))
