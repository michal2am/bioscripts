import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
from scipy.interpolate import interp1d
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
concentrations = [0.001, 0.001, 0.01, 0.1, 1, 10, 100]
models = []

all_time = tm.time()


for conc in concentrations:

    '''
    stimulus = Stimulus.pair_square([0., 10.], [5., 15.], 10.0)
    suspend = [[5., 10.], [15., 20.]]
    '''

    stimulus = Stimulus.square(0., 25., conc)
    suspend = [[2.5], [10.]]


    log.info("### Adding rates:")

    if model == 'kisiel':

        r_kon   = Kinetic.State.Rate('kon',     2*1e7*1e-6,   stimulus)
        r_2kon  = Kinetic.State.Rate('2kon',    2*2*1e7*1e-6, stimulus)
        r_koff  = Kinetic.State.Rate('koff',    1116*1e-3)
        r_2koff = Kinetic.State.Rate('2koff',   2*1116*1e-3)
        r_b1    = Kinetic.State.Rate('b1',      150*1e-3)
        r_b2    = Kinetic.State.Rate('b2',      18000*1e-3)
        r_b3    = Kinetic.State.Rate('b3',      35*1e-3)
        r_b4    = Kinetic.State.Rate('b4',      200*1e-3)
        r_b5    = Kinetic.State.Rate('b5',      100*1e-3)
        r_a1    = Kinetic.State.Rate('a1',      20000*1e-3)
        r_a2    = Kinetic.State.Rate('a2',      800*1e-3)
        r_a3    = Kinetic.State.Rate('a3',      3333*1e-3)
        r_a4    = Kinetic.State.Rate('a4',      17500*1e-3)
        r_a5    = Kinetic.State.Rate('a5',      3000*1e-3)
        r_d1    = Kinetic.State.Rate('d1',      310*1e-3)
        r_d2    = Kinetic.State.Rate('d2',      800*1e-3)
        r_d4    = Kinetic.State.Rate('d4',      1000*1e-3)
        r_r1    = Kinetic.State.Rate('r1',      5*1e-3)
        r_r2    = Kinetic.State.Rate('r2',      400*1e-3)
        r_r4    = Kinetic.State.Rate('r4',      10*1e-3)
        r_y2    = Kinetic.State.Rate('y2',      4400*1e-3)
        r_g2    = Kinetic.State.Rate('g2',      4300*1e-3)
        r_0     = Kinetic.State.Rate('block',   0.00*1e-3)

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
        st_a1o  = Kinetic.State(0, 'A1O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_a1,   r_0,    r_0,    r_d1,   r_0,    r_0],   0, 'open')
        st_a2o  = Kinetic.State(1, 'A2O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_a2,   r_0,    r_0,    r_0],   0, 'open')
        st_a3o  = Kinetic.State(2, 'A3O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_a3,   r_0,    r_0],   0, 'open')
        st_a4o  = Kinetic.State(3, 'A4O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_a4,   r_0,    r_0,    r_0,    r_d4],  0, 'open')
        st_a5o  = Kinetic.State(4, 'A5O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_a5],  0, 'open')
        st_r    = Kinetic.State(5, 'R'  , [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_2kon, r_0,    r_0,    r_0,    r_0,    r_0],   1, 'unbound')
        st_a1r  = Kinetic.State(6, 'A1R', [r_a1,    r_0,    r_0,    r_0,    r_0,    r_koff, r_0,    r_kon,  r_0,    r_0,    r_0,    r_0],   0, 'single-bound')
        st_a2r  = Kinetic.State(7, 'A2R', [r_0,     r_0,    r_0,    r_b4,   r_0,    r_0,    r_2koff,r_0,    r_g2,   r_0,    r_0,    r_0],   0, 'double-bound')
        st_a2f  = Kinetic.State(8, 'A2F', [r_0,     r_b2,   r_0,    r_0,    r_0,    r_0,    r_0,    r_y2,    r_0,   r_0,    r_d2,   r_0],   0, 'flipped')
        st_a1d  = Kinetic.State(9, 'A1D', [r_r1,    r_0,    r_b3,   r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0],   0, 'desensitized')
        st_a2d  = Kinetic.State(10,'A2D', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_r2,   r_0,    r_0,    r_0],   0, 'desensitized')
        st_a4d  = Kinetic.State(11,'A4D', [r_0,     r_0,    r_0,    r_r4,   r_b5,   r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0],   0, 'desensitized')

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

    models.append(model_kinetic)


# solving!

solve_ode = False
solve_ode_dr = True
solve_glp = False
solve_sto = False
plot_stim = True


part = 10
t0 = 0
te = 100

# figure preparation

sns.set_style("ticks", {'legend.frameon': True})
sns.set_context("poster")
sns.set_palette('Paired', states_no)
colors = sns.color_palette('deep', states_no)

if solve_ode:

    fig_ode = plt.figure()
    fig_ode.canvas.set_window_title('ODE results')
    gs = grd.GridSpec(6, 2, width_ratios=[3, 2], hspace=2, wspace=0.20)
    ax2 = plt.subplot(gs[5, 0])
    ax1 = plt.subplot(gs[0:5, 0], sharex=ax2)
    ax4 = plt.subplot(gs[3:, 1])
    ax3 = plt.subplot(gs[0:3, 1])
    fig_ode.add_subplot(ax1, ax2, ax3, ax4)

    # ax1 trajectory
    solver = SolverOde(model_kinetic.trmn, ini_conc, states_na, t0, te, False)
    TP = solver.get_results()
    TP.plot(x='time', ax=ax1, legend=True)

    open_states = [state_model.name for state_model in states if state_model.category == 'open']
    cumulative_open = TP[open_states].sum(axis=1)
    cumulative_open = pd.concat([TP['time'], cumulative_open], axis=1)
    cumulative_open.columns = ['time', 'AxO']
    cumulative_open.plot(x='time', ax=ax1, legend=True, color='k')

    ax1.legend(loc='upper center', ncol=5, borderaxespad=0.)
    ax1.set_xlim([t0 - 0.5, te + 0.5])
    ax1.set_ylim([0, 1.00])
    ax1.set_ylabel('state probability')

    # ax2 stimulus
    T = np.linspace(t0 - 2 , te, 1e5)
    stim = [stimulus(ti) for ti in T]
    ax2.step(T, stim, 'k-')

    ax2.set_xlim([t0 - 5, te + 2])
    ax2.set_ylim([0, max(stim)])
    ax2.set_yticks([0, max(stim)])

    ax2.set_ylabel('stimulus [mM]')
    ax2.set_xlabel('time [ms]')

    # ax3 equilibrium
    solver = SolverOde(model_kinetic.trmn, ini_conc, states_na, t0, te, True)
    TP = solver.get_results()[-1:].drop('time', 1)
    log.info('###Equilibrium by states:')
    log.info(TP)
    sns.barplot(data=TP, ax=ax3)

    ax3.set_ylabel('equilibrium occupancies')
    ax3.set_xlabel('state')
    ax3.set(yscale='log')

    # ax4 equilibrium categorical
    categories = pd.DataFrame({'open': [0], 'unbound': [0], 'single-bound': [0], 'double-bound': [0], 'flipped': [0], 'desensitized': [0]})
    for category in categories.columns:
        for state_model in states:
            if state_model.category == category:
                categories[state_model.category].iloc[0] += TP[state_model.name].iloc[0]

    log.info('###Equilibrium by categories:')
    log.info(categories)
    categories = categories.transpose()
    categories.columns = ['cumulative equilibrium occupancy']
    categories.plot.pie(y='cumulative equilibrium occupancy', legend=True, labels=None, ax=ax4)

    ax4.legend(loc='lower center', nrow=1)
    ax4.axis('equal')

if solve_ode_dr:

    fig_ode = plt.figure()
    fig_ode.canvas.set_window_title('ODE results')
    gs = grd.GridSpec(2, 2, width_ratios=[1, 1], hspace=.1, wspace=0.)
    ax2 = plt.subplot(gs[1, 0:])
    ax1 = plt.subplot(gs[0, 0:], sharex=ax2)


    fig_ode.add_subplot(ax1)

    concs = pd.DataFrame({'concentration': concentrations})

    solver = SolverOde(models[0].trmn, ini_conc, states_na, t0, te, True)

    # first by state
    cumulative = solver.get_results()[-1:].drop('time', 1)

    # first by category
    categories = pd.DataFrame({'open': [0], 'unbound': [0], 'single-bound': [0], 'double-bound': [0], 'flipped': [0], 'desensitized': [0]})
    for category in categories.columns:
        for state_model in states:
            if state_model.category == category:
                categories[state_model.category].iloc[0] += cumulative[state_model.name].iloc[0]


    for model in models[1:]:
        solver = SolverOde(model.trmn, ini_conc, states_na, t0, te, True)
        TP = solver.get_results()[-1:].drop('time', 1)
        log.info(TP)
        cumulative = cumulative.append(TP)

        new_categories = pd.DataFrame({'open': [0], 'unbound': [0], 'single-bound': [0], 'double-bound': [0], 'flipped': [0], 'desensitized': [0]})
        for category in categories.columns:
            for state_model in states:
                if state_model.category == category:
                    new_categories[state_model.category].iloc[0] += TP[state_model.name].iloc[0]
        categories = categories.append(new_categories)

    cumulative = cumulative.reset_index(drop=True)
    cumulative = pd.concat([concs, cumulative], axis=1)
    categories = categories.reset_index(drop=True)
    categories = pd.concat([concs, categories], axis=1)

    #conc_fit = np.linspace(concentrations[0], concentrations[-1], num=100)
    #resp_fit = interp1d(categories.as_matrix(columns='concentration'), categories.as_matrix(columns='open'), kind='cubic')
    #dr_fit = resp_fit(conc_fit)
    #print(dr_fit)

    cumulative.plot(x='concentration', ax=ax1, style='o')
    categories.plot(x='concentration', ax=ax2, style='o')
    #ax2.plot(conc_fit, dr_fit)

    log.info('###Equilibrium by states:')
    log.info(cumulative)
    log.info('###Equilibrium by categories:')
    log.info(categories)

    ax1.set(xscale='log')
    ax1.set_ylabel("equilibrium occupancy")
    ax2.set(xscale='log')
    ax2.set_xlabel("stimulus [mM]")
    ax2.set_ylabel("equilibrium occupancy")
    ax1.legend(loc='upper center', ncol=6, borderaxespad=0.)
    ax2.legend(loc='upper center', ncol=6, borderaxespad=0.)
    ax2.set_xlim(-10, 110)

    # ax3 distributions
    # lifetimes = pd.DataFrame()
    # time = np.linspace(0., 100, num=1000)
    # lifetimes['time'] = time
    # sumrates = [np.sum(state_rates) for state_rates in model_kinetic.trm(0)]
    # for no, name in enumerate(states_na):
    #     lifetimes[name] = sts.expon.pdf(time, loc=0, scale=1/sumrates[no])
    # lifetimes.plot(x='time', ax=ax3, legend=False, logx=True, logy=False)
    # ax3.set_ylabel('probability density function')
    # ax3.set_xlabel('occupancy time [ms]')



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
