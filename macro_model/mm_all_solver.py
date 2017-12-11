import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.stats as sts
import time as tm
import logging as log
from mm_solver_glp import SolverGlp
from mm_kinetic_models import ModelBuilder
from mm_analyzer_ode import AnalyzerODE

log.basicConfig(filename='mm_kin.log', filemode='w', format='%(message)s', level=log.DEBUG)
log.info("### Kinetic solver starts!")
# log.captureWarnings(True)
all_time = tm.time()

# model configuration
build_models = ModelBuilder('fjwm', [0.001, 0.01, 0.1, 1, 10], 'single')
t0 = 0
te = 100
part = 10

# solving!
solve_ode = True
dynamic = True
steady = True

solve_glp = False

if solve_ode:
    ode_analysis = AnalyzerODE(build_models, t0, te)
    if dynamic:
        ode_analysis.plot_dynamic_response(3)
    if steady:
        ode_analysis.plot_steady_dose_response()

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


plt.show()
print("--- %s seconds ---" % (tm.time() - all_time))
