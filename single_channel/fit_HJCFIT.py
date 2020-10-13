import sys
import time
import math
import numpy as np
import pandas as pd
import argparse
from scipy.optimize import minimize

from dcpyps.samples import samples
from dcpyps import dataset
from dcpyps import mechanism

from dcprogs.likelihood import Log10Likelihood

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
parser.add_argument('--files', nargs='+')
parser.add_argument('--tres', type=float)
parser.add_argument('--tcrit', type=float)
parser.add_argument('--config')
parser.add_argument('--project')
args = parser.parse_args()

if args.config:
    config = pd.read_csv(args.config)
    results = []

    for file_name in config.file.unique():

        single_cell = config[config.loc[:, 'file'] == file_name].copy()
        single_cell.reset_index(inplace=True)

        sc_model = single_cell.at[0, 'model']
        sc_type = single_cell.at[0, 'type']
        sc_tres = single_cell.at[0, 'tres']/1000000
        sc_tcrit = single_cell.at[0, 'tcrit']/1000
        sc_scns = list(single_cell.loc[:, 'file_scn'])
        sc_scns = [name + '.SCN' for name in sc_scns]
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
        likelihood = Log10Likelihood(bursts, mec.kA, sc_tres, sc_tcrit)

        iternum = 0
        start = time.clock()
        res = minimize(dcprogslik, np.log(theta), method='Nelder-Mead', callback=printiter, )
        t3 = time.clock() - start
        print("\n\n\nScyPy.minimize (Nelder-Mead) Fitting finished: %4d/%02d/%02d %02d:%02d:%02d\n"
              % time.localtime()[0:6])
        print('time in ScyPy.minimize (Nelder-Mead)=', t3)
        print('xout', res.x)
        mec.theta_unsqueeze(np.exp(res.x))
        print("\n Final rate constants:")
        mec.printout(sys.stdout)
        lik = dcprogslik(res.x)
        print("\nFinal likelihood = {0:.6f}".format(-lik))

        if sc_model == 'CO':

            alpha = mec.Rates[0].rateconstants[0]
            beta = mec.Rates[1].rateconstants[0]

            refer_result = {'project': args.project, 'type': sc_type, 'file': file_name, 'model': sc_model,
                            'equilibrium': np.log(beta / alpha), 'forward': np.log(beta),
                            'alpha': alpha, 'beta': beta,
                            }

        elif sc_model == 'CFO':

            alpha = mec.Rates[0].rateconstants[0]
            beta = mec.Rates[1].rateconstants[0]
            gamma = mec.Rates[2].rateconstants[0]
            delta = mec.Rates[3].rateconstants[0]

            refer_result = {'project': args.project, 'type': sc_type, 'file': file_name, 'model': sc_model,
                            'equilibrium': np.log((delta * beta) / (gamma * alpha)),
                            'forward': np.log(delta * beta / gamma),
                            'alpha': alpha, 'beta': beta,
                            'delta': delta, 'gamma': gamma,
                            }


        results.append(refer_result)

    results = pd.DataFrame(results)
    print(results)
    results.to_csv(args.project + '_hjcfit_rates_all.csv')


else:

    ### running without config file, rather legacy for testing only, probably refecator/update needed

    rec = dataset.SCRecord(args.files, 100e-9, args.tres, args.tcrit)
    rec.record_type = 'recorded'
    rec.printout()

    mec = mechanism_RO([5000, 5000])

    mec.update_mr()
    mec.printout(sys.stdout)
    theta = mec.theta()
    print ('\ntheta=', theta)

    bursts = rec.bursts.intervals()
    logfac = math.log(10)
    likelihood = Log10Likelihood(bursts, mec.kA, args.tres, args.tcrit)

    iternum = 0
    start = time.clock()
    res = minimize(dcprogslik, np.log(theta), method='Nelder-Mead', callback=printiter,)
    t3 = time.clock() - start
    print ("\n\n\nScyPy.minimize (Nelder-Mead) Fitting finished: %4d/%02d/%02d %02d:%02d:%02d\n"
        %time.localtime()[0:6])
    print ('time in ScyPy.minimize (Nelder-Mead)=', t3)
    print ('xout', res.x)
    mec.theta_unsqueeze(np.exp(res.x))
    print ("\n Final rate constants:")
    mec.printout(sys.stdout)
    lik = dcprogslik(res.x)
    print ("\nFinal likelihood = {0:.6f}".format(-lik))

