import sys
import time
import math
import numpy as np
import pandas as pd
import argparse
import re
from scipy.optimize import minimize

from dcpyps.samples import samples
from dcpyps import dataset
from dcpyps import mechanism

from dcprogs.likelihood import Log10Likelihood
from scalcs import scalcslib as scl
from scalcs import scplotlib as scpl
import matplotlib.pyplot as plt




def dcprogslik(x):
    mec.theta_unsqueeze(np.exp(x))
    mec.set_eff('c', 100e-9)
    return -likelihood(mec.Q) * logfac


def printiter(theta):
    global iternum
    iternum += 1
    lik = dcprogslik(theta)
    if iternum % 10 == 0:
        print("iteration # {0:d}; log-lik = {1:.6f}".format(iternum, -lik))
        print(np.exp(theta))


def mechanism_RO(rates):

    o = mechanism.State('A', 'F*', 50e-12)
    r = mechanism.State('B', 'R', 0.0)

    rate_list = [
        mechanism.Rate(rates[1], o, r, name='alpha', limits=[1e-15, 1e+7]),
        mechanism.Rate(rates[0], r, o, name='beta', limits=[1e-15, 1e+7]),
    ]

    complete_mechanism = mechanism.Mechanism(rate_list, mtitle='CO', rtitle='CO_rates')
    complete_mechanism.set_eff('c', 100e-9)

    return complete_mechanism


def mechanism_RFO(rates):

    o = mechanism.State('A', 'F*', 50e-12)
    f = mechanism.State('B', 'F', 0.0)
    r = mechanism.State('C', 'R', 0.0)

    rate_list = [
        mechanism.Rate(rates[3], o, f, name='alpha', limits=[1e+0, 1.5e+4]),
        mechanism.Rate(rates[2], f, o, name='beta', limits=[1e+0, 1.5e+4]),
        mechanism.Rate(rates[1], f, r, name='gamma', limits=[1e+0, 1.5e+4]),
        mechanism.Rate(rates[0], r, f, name='delta', limits=[1e+0, 1.5e+4])
    ]

    complete_mechanism = mechanism.Mechanism(rate_list, mtitle='CFO', rtitle='CFO_rates')
    complete_mechanism.set_eff('c', 100e-9)

    return complete_mechanism


parser = argparse.ArgumentParser()
parser.add_argument('--config')
args = parser.parse_args()

config = pd.read_csv(args.config)
project = '_'.join(re.split('[_,.]', args.config)[2:-1])
results = []

for file_name in config.file.unique():

    single_cell = config[config.loc[:, 'file'] == file_name].copy()
    single_cell.reset_index(inplace=True)

    sc_model = single_cell.at[0, 'model']
    sc_type = single_cell.at[0, 'type']
    sc_tres = single_cell.at[0, 'tres']/1000000
    sc_tcrit = single_cell.at[0, 'tcrit']/1000

    # TODO: temporary commented out for CO testings
    # sc_t1_exp = single_cell.at[0, 't1_exp']
    # sc_t2_exp = single_cell.at[0, 't2_exp']
    # sc_p1_exp = single_cell.at[0, 'p1_exp']
    #sc_p2_exp = single_cell.at[0, 'p2_exp']

    sc_scns = list(single_cell.loc[:, 'file_scn'])
    sc_scns = [name if name.endswith('.SCN') else name + '.SCN' for name in sc_scns]
    # TODO: should upper only the .scn extension
    sc_scns = [full_name.upper() for full_name in sc_scns]
    rec = dataset.SCRecord(sc_scns, 100e-9, sc_tres, sc_tcrit)
    rec.record_type = 'recorded'
    rec.printout()

    if sc_model == 'CO':
        mec = mechanism_RO([single_cell.at[0, 'beta'], single_cell.at[0, 'alpha']])
    elif sc_model == 'CFO':
        mec = mechanism_RFO([single_cell.at[0, 'delta'], single_cell.at[0, 'gamma'],
                             single_cell.at[0, 'beta'], single_cell.at[0, 'alpha']])

    mec.update_mr()
    mec.printout(sys.stdout)
    theta = mec.theta()
    print('\ntheta=', theta)

    bursts = rec.bursts.intervals()
    logfac = math.log(10)
    print('\nFirst likelihood calculation ...')
    print(sc_tres, sc_tcrit)
    likelihood = Log10Likelihood(bursts, mec.kA, sc_tres, sc_tcrit)

    print('\nFirst iter ...')
    iternum = 0
    #start = time.clock()
    print(dcprogslik(theta))
    print("Minimizing starts")
    res = minimize(dcprogslik, np.log(theta), method='Powell', callback=printiter, )
    #t3 = time.clock() - start
    #print("\n\n\nScyPy.minimize (Nelder-Mead) Fitting finished: %4d/%02d/%02d %02d:%02d:%02d\n"
    #      % time.localtime()[0:6])
    #print('time in ScyPy.minimize (Nelder-Mead)=', t3)
    print('xout', res.x)
    mec.theta_unsqueeze(np.exp(res.x))
    print("\n Final rate constants:")
    mec.printout(sys.stdout)
    lik = dcprogslik(res.x)
    print("\nFinal likelihood = {0:.6f}".format(-lik))

    # TODO: temporary commented out for CO testings

    '''
    print(scl.printout_distributions(mec, sc_tres))
    shuts = scl.printout_distributions(mec, sc_tres)

    t1_mod = shuts.splitlines()[10].split('\t')[1]
    p1_mod = shuts.splitlines()[10].split('\t')[2]
    t2_mod = shuts.splitlines()[11].split('\t')[1]
    p2_mod = shuts.splitlines()[11].split('\t')[2]

    shuts_format = 't1: ' + t1_mod + ' p1: ' + p1_mod + ' t2: ' + t2_mod + ' p2: ' + p2_mod

    t, ipdf, epdf, apdf = scpl.shut_time_pdf(mec, sc_tres)
    plt.semilogx(t, ipdf, 'r--', t, epdf, 'b-', t, apdf, 'g-')
    plt.ylabel('fshut(t)')
    plt.xlabel('Shut time, ms')
    plt.title(file_name + ' ' + sc_type + ' ' + str(sc_tres*1000000) + ' ' + str(sc_tcrit*1000) + '\n' + shuts_format)
    print('RED- ideal distribution\nGREEN- HJC distribution (corrected for missed events)')
    plt.savefig(project + '_' + file_name.strip('.abf') + '_shut_plot.png')
    plt.close()
    #plt.show()
    '''

    if sc_model == 'CO':

        alpha = mec.Rates[0].rateconstants[0]
        beta = mec.Rates[1].rateconstants[0]

        refer_result = {'project': project, 'type': sc_type, 'file': file_name, 'model': sc_model,
                        'equilibrium': np.log(beta / alpha), 'forward': np.log(beta),
                        'alpha': alpha, 'beta': beta,
                        }

    elif sc_model == 'CFO':

        alpha = mec.Rates[0].rateconstants[0]
        beta = mec.Rates[1].rateconstants[0]
        gamma = mec.Rates[2].rateconstants[0]
        delta = mec.Rates[3].rateconstants[0]

        refer_result = {'project': project, 'type': sc_type, 'file': file_name, 'model': sc_model,
                        # 'f1': np.log(delta),
                        # 'f2': np.log(delta/gamma),
                        # 'f3': np.log(delta * beta / gamma),
                        # 'e': np.log((delta * beta) / (gamma * alpha)),
                        'alpha': alpha, 'beta': beta,
                        'gamma': gamma, 'delta': delta,
                        't1_mod': t1_mod, 'p1_mod': p1_mod, 't2_mod': t2_mod, 'p2_mod': p2_mod,
                        't1_exp': sc_t1_exp, 'p1_exp': sc_p1_exp, 't2_exp': sc_t2_exp, 'p2_exp': sc_p2_exp
                        }


    results.append(refer_result)

results = pd.DataFrame(results)
print(results)
results.to_csv('hjcfit_rates_' + project + '.csv')


