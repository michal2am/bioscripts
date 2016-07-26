import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import matplotlib.colors as cls
from collections import OrderedDict
from scipy.interpolate import interp1d
import seaborn as sns
import pandas as pd
import scipy.stats as sts
import time as tm
import logging as log
from mm_solver_ode import SolverOde
from mm_solver_glp import SolverGlp
from mm_kinetic_models import ModelBuilder

log.basicConfig(filename='mm_kin.log', filemode='w', format='%(message)s', level=log.DEBUG)
log.info("### Kinetic solver starts!")


all_time = tm.time()

build_models = ModelBuilder('kisiel', [0.001, 0.01, 0.1, 1, 10, 100], 'single')
models = build_models.models
agonist_concentrations = build_models.agonist_concentrations

states = build_models.models[0].states
states_names = build_models.models[0].states_names
states_categories = build_models.models[0].states_categories
states_belongs = build_models.models[0].states_belongs
states_number = build_models.models[0].states_number
states_ini_concentrations = build_models.models[0].states_ini_concentrations


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
sns.set_palette('Paired', states_number)
colors = sns.color_palette('Paired', states_number)
for color in colors:
    print(color)

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

    fig_ode_dr = plt.figure()
    fig_ode_dr.canvas.set_window_title('ODE dose response results')
    gs = grd.GridSpec(2, 2, width_ratios=[1, 1], hspace=.1, wspace=0.)
    ax2 = plt.subplot(gs[1, 0:])
    ax1 = plt.subplot(gs[0, 0:], sharex=ax2)
    fig_ode_dr.add_subplot(ax1, ax2)

    ax1.set(xscale='log')
    ax1.set_ylabel("equilibrium occupancy")
    ax2.set(xscale='log')
    ax2.set_xlabel("stimulus [mM]")
    ax2.set_ylabel("equilibrium occupancy")

    ax2.set_xlim(-10, 110)

    concs = pd.DataFrame({'concentration': agonist_concentrations})
    dr_bystate = pd.DataFrame()

    for model in models:
        solver = SolverOde(model.trmn, states_ini_concentrations, states_names, t0, te, True)
        dr_bystate = dr_bystate.append(solver.get_results()[-1:].drop('time', 1))

    dr_bystate = dr_bystate.reset_index(drop=True)
    dr_bystate = pd.concat([concs, dr_bystate], axis=1)
    dr_bystate = dr_bystate.set_index('concentration', drop=True)

    dr_bycategory = pd.DataFrame({category: dr_bystate[states_belongs[category]].sum(axis=1) for category in states_belongs})

    dr_bystate_fit = dr_bystate.copy()
    dr_bycategory_fit = dr_bycategory.copy()

    def fit_logarithmic(data):
        """
        fits numeric interpolation with logarithmic x-values conversion
        :param data: pandas dataframe, x-values as index
        :return: pandas dataframe, interpolated x-values as index
        """

        xlabel = data.index.name
        xknown = data.index

        # replace index with index logarithm
        data = data.reset_index()
        data[xlabel] = data[xlabel].apply(np.log10)
        data = data.set_index(xlabel, drop=True)

        # prepare extended index to interpolate
        log_xknown = np.log10(xknown)
        log_xunknown = np.linspace(log_xknown[0], log_xknown[-1], num=100)
        log_x = np.concatenate((log_xknown[1:-1], log_xunknown), axis=0)
        log_x = np.sort(log_x)

        # insert new index and interpolate
        data = data.reindex(log_x)
        data = data.interpolate(method='akima')

        # invert index back
        data = data.reset_index()
        data[xlabel] = data[xlabel].rpow(np.ones(len(data[xlabel]))*10)
        data = data.set_index(xlabel, drop=True)

        return data

    dr_bystate_fit = fit_logarithmic(dr_bystate_fit)
    dr_bycategory_fit = fit_logarithmic(dr_bycategory_fit)

    dr_bystate.plot(ax=ax1, marker='o', linestyle='None')
    dr_bystate_fit.plot(ax=ax1, legend=False)
    dr_bycategory.plot(ax=ax2, style='o', colormap=cls.ListedColormap(colors[0:6]))
    dr_bycategory_fit.plot(ax=ax2, legend=False, colormap=cls.ListedColormap(colors[0:6]))

    ax1.legend(numpoints=1)

    #ax1.legend(loc='upper center', ncol=6, borderaxespad=0.)
    #ax2.legend(loc='upper center', ncol=6, borderaxespad=0.)

    log.info('###Equilibrium by states:')
    log.info(dr_bystate)
    log.info('###Equilibrium by categories:')
    log.info(dr_bycategory)



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
