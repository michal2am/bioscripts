import sys
import math
import numpy as np
import pandas as pd
import argparse
import re
import os
import matplotlib.pyplot as plt

from scipy.optimize import minimize
from dcpyps import dataset
from dcpyps import mechanism
from dcprogs.likelihood import Log10Likelihood
from scalcs import scalcslib as scl
from scalcs import scplotlib as scpl



def dcprogslik(x):
    mec.theta_unsqueeze(np.exp(x))
    mec.set_eff('c', 100e-9)
    return -likelihood(mec.Q) * logfac

def printiter(theta):
    global iternum
    iternum += 1
    lik = dcprogslik(theta)
    if iternum % 10 == 0:
        pass
        #print("iteration # {0:d}; log-lik = {1:.6f}".format(iternum, -lik))
        #print(np.exp(theta))

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

def mechanism_CFOODD(rates):

    c = mechanism.State('B', 'A2R', 0.0)
    f = mechanism.State('B', 'A2F', 0.0)
    o = mechanism.State('A', 'A2O', 50e-12)
    op = mechanism.State('A', 'A2Op', 50e-12)
    d = mechanism.State('C', 'A2D', 0.0)
    dp = mechanism.State('C', 'A2Dp', 0.0)

    rate_list = [
        mechanism.Rate(rates[0], f, o, name='beta', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[1], f, op, name='betap', limits=[1e0, 1.5e4]), #1e0
        mechanism.Rate(rates[2], o, f, name='alpha', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[3], op, f, name='alphap', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[4], c, f, name='delta', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[5], f, c, name='gamma', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[6], f, d, name='d', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[7], f, dp, name='dp', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[8], d, f, name='r', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[9], dp, f, name='rp', limits=[1e0, 1.5e4])
    ]

    complete_mechanism = mechanism.Mechanism(rate_list, mtitle='CFOODD', rtitle='CFOODD_rates')
    #complete_mechanism.Rates[1].fixed = True
    # complete_mechanism.set_eff('c', 100e-9)

    return complete_mechanism
def mechanism_CFODD(rates):

    c = mechanism.State('B', 'A2R', 0.0)
    f = mechanism.State('B', 'A2F', 0.0)
    o = mechanism.State('A', 'A2O', 50e-12)
    d = mechanism.State('C', 'A2D', 0.0)
    dp = mechanism.State('C', 'A2Dp', 0.0)

    rate_list = [
        mechanism.Rate(rates[0], f, o, name='beta', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[1], o, f, name='alpha', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[2], c, f, name='delta', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[3], f, c, name='gamma', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[4], f, d, name='d', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[5], f, dp, name='dp', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[6], d, f, name='r', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[7], dp, f, name='rp', limits=[1e0, 1.5e4])
    ]

    complete_mechanism = mechanism.Mechanism(rate_list, mtitle='CFOODD', rtitle='CFOODD_rates')
    #complete_mechanism.Rates[1].fixed = True
    # complete_mechanism.set_eff('c', 100e-9)

    return complete_mechanism

def mechanism_CFOOD(rates):

    c = mechanism.State('B', 'A2R', 0.0)
    f = mechanism.State('B', 'A2F', 0.0)
    o = mechanism.State('A', 'A2O', 50e-12)
    op = mechanism.State('A', 'A2Op', 50e-12)
    d = mechanism.State('C', 'A2D', 0.0)

    rate_list = [
        mechanism.Rate(rates[0], f, o, name='beta', limits=[1e0, 1e1]),
        mechanism.Rate(rates[1], f, op, name='betap', limits=[1e0, 1e1]),
        mechanism.Rate(rates[2], o, f, name='alpha', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[3], op, f, name='alphap', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[4], c, f, name='delta', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[5], f, c, name='gamma', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[6], f, d, name='d', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[7], d, f, name='r', limits=[1e0, 1.5e4]),
    ]

    complete_mechanism = mechanism.Mechanism(rate_list, mtitle='CFOODD', rtitle='CFOODD_rates')
    # complete_mechanism.set_eff('c', 100e-9)

    return complete_mechanism

def mechanism_CFOD(rates):

    c = mechanism.State('B', 'A2R', 0.0)
    f = mechanism.State('B', 'A2F', 0.0)
    o = mechanism.State('A', 'A2O', 50e-12)
    d = mechanism.State('C', 'A2D', 0.0)

    rate_list = [
        mechanism.Rate(rates[0], f, o, name='beta', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[1], o, f, name='alpha', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[2], c, f, name='delta', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[3], f, c, name='gamma', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[4], f, d, name='d', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[5], d, f, name='r', limits=[1e0, 1.5e4]),
    ]

    complete_mechanism = mechanism.Mechanism(rate_list, mtitle='CFOODD', rtitle='CFOODD_rates')
    # complete_mechanism.set_eff('c', 100e-9)

    return complete_mechanism

def mechanism_CFOO(rates):

    c = mechanism.State('B', 'A2R', 0.0)
    f = mechanism.State('B', 'A2F', 0.0)
    o = mechanism.State('A', 'A2O', 50e-12)
    op = mechanism.State('A', 'A2Op', 50e-12)

    rate_list = [
        mechanism.Rate(rates[0], f, o, name='beta', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[1], f, op, name='betap', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[2], o, f, name='alpha', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[3], op, f, name='alphap', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[4], c, f, name='delta', limits=[1e0, 1.5e4]),
        mechanism.Rate(rates[5], f, c, name='gamma', limits=[1e0, 1.5e4]),
    ]

    complete_mechanism = mechanism.Mechanism(rate_list, mtitle='CFOODD', rtitle='CFOODD_rates')
    # complete_mechanism.set_eff('c', 100e-9)

    return complete_mechanism

parser = argparse.ArgumentParser()
parser.add_argument('--config')
parser.add_argument('--tcrit')
parser.add_argument('--project')
args = parser.parse_args()

config = pd.read_csv(args.config)
#project = '_'.join(re.split('[_,.]', args.config)[2:-1])
project = args.project
results = []
model_results = []

for file_name in config.file.unique():
    # TODO: BIG PROBLEM IF MULTIPLE ABFs FOR SINGLE CELL

    single_cell = config[config.loc[:, 'file'] == file_name].copy()
    single_cell.reset_index(inplace=True)

    sc_model = single_cell.at[0, 'model']
    sc_type = single_cell.at[0, 'type']
    sc_tres = (single_cell.at[0, 'tres']/1000000)
    sc_tcrit = (single_cell.at[0, args.tcrit]/1000)


    sc_scns = list(single_cell.loc[:, 'file_scn'])
    sc_scns = [name if name.endswith('.SCN') else name + '.SCN' for name in sc_scns]
    print(sc_scns)

    checked_scns = []
    for scn in sc_scns:
        if scn in (os.listdir()):
            print('SCN file found {}.'.format(scn))
            checked_scns.append(scn)
        elif scn.upper() in (os.listdir()):
            print('SCN file found {}, changing to uppercase.'.format(scn))
            up_scn = scn.upper()
            checked_scns.append(up_scn)
        else:
            print('SCN file not found {}.'.format(scn))

    # TODO: tcrit for non-CO models
    if sc_model == 'CFOODD':
        rec = dataset.SCRecord(checked_scns, 100e-9, sc_tres, 1000000)
    elif sc_model == 'CFODD':
        rec = dataset.SCRecord(checked_scns, 100e-9, sc_tres, 1000000)
    else:
        rec = dataset.SCRecord(checked_scns, 100e-9, sc_tres, sc_tcrit)

    rec.record_type = 'recorded'
    rec.printout()

    if sc_model == 'CO':
        mec = mechanism_RO([single_cell.at[0, 'beta'], single_cell.at[0, 'alpha']])
    elif sc_model == 'CFO':
        mec = mechanism_RFO([single_cell.at[0, 'delta'], single_cell.at[0, 'gamma'],
                             single_cell.at[0, 'beta'], single_cell.at[0, 'alpha']])
    elif sc_model == 'CFOO':
        mec = mechanism_CFOO([single_cell.at[0, 'beta'], single_cell.at[0, 'betap'],
                          single_cell.at[0, 'alpha'], single_cell.at[0, 'alphap'],
                          single_cell.at[0, 'delta'], single_cell.at[0, 'gamma']])
    elif sc_model == 'CFOOD':
        mec = mechanism_CFOOD([single_cell.at[0, 'beta'], single_cell.at[0, 'betap'],
                          single_cell.at[0, 'alpha'], single_cell.at[0, 'alphap'],
                          single_cell.at[0, 'delta'], single_cell.at[0, 'gamma'],
                          single_cell.at[0, 'd'],
                          single_cell.at[0, 'r']])
    elif sc_model == 'CFOD':
        mec = mechanism_CFOOD([single_cell.at[0, 'beta'],
                          single_cell.at[0, 'alpha'],
                          single_cell.at[0, 'delta'], single_cell.at[0, 'gamma'],
                          single_cell.at[0, 'd'],
                          single_cell.at[0, 'r']])
    elif sc_model == 'CFODD':
        print('making mechanism')
        mec = mechanism_CFODD([single_cell.at[0, 'beta'],
                          single_cell.at[0, 'alpha'],
                          single_cell.at[0, 'delta'], single_cell.at[0, 'gamma'],
                          single_cell.at[0, 'd'], single_cell.at[0, 'dp'],
                          single_cell.at[0, 'r'], single_cell.at[0, 'rp']])
    elif sc_model == 'CFOODD':
        print('making mechanism')
        mec = mechanism_CFOODD([single_cell.at[0, 'beta'], single_cell.at[0, 'betap'],
                          single_cell.at[0, 'alpha'], single_cell.at[0, 'alphap'],
                          single_cell.at[0, 'delta'], single_cell.at[0, 'gamma'],
                          single_cell.at[0, 'd'], single_cell.at[0, 'dp'],
                          single_cell.at[0, 'r'], single_cell.at[0, 'rp']])

    mec.update_mr()
    mec.printout(sys.stdout)
    theta = mec.theta()

    bursts = rec.bursts.intervals()
    logfac = math.log(10)
    print('\nFirst likelihood calculation ...')
    print(sc_tres, sc_tcrit)
    likelihood = Log10Likelihood(bursts, mec.kA, sc_tres, sc_tcrit)

    print('\nFirst iter ...')
    iternum = 0

    print("Minimizing starts")
    res = minimize(dcprogslik, np.log(theta), method='Nelder-Mead', callback=printiter, )

    print('xout', res.x)
    mec.theta_unsqueeze(np.exp(res.x))
    print("\n Final rate constants:")
    mec.printout(sys.stdout)
    lik = dcprogslik(res.x)
    print("\nFinal likelihood = {0:.6f}".format(-lik))


    plots = True
    if plots:

        # plot event times of shuts, works only for CFO model
        if sc_model == 'CFO':
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


        # plot and events times for openings and shuts, works only for complex models
        if sc_model in ['CFOODD', 'CFOOD', 'CFOO', 'CFODD']:

            (open_dist, asymptotic_open_dist, asymptotic_open_norm,
             shut_dist, asymptotic_shut_dist, asymptotic_shut_norm,
             event_times_data) = scl.printout_distributions(mec, sc_tres)

            print(event_times_data)

            open_dist = pd.DataFrame(open_dist, columns=['term', 'w', 'rate', 'op_tau_id', 'op_area_id']).drop(['w', 'rate'], axis=1)
            asymptotic_open_dist = pd.DataFrame(asymptotic_open_dist, columns=['term', 'op_tau_as', 'op_area_as', 'rate']).drop(['rate'], axis=1)
            asymptotic_open_norm = pd.DataFrame(asymptotic_open_norm, columns=['term', 'op_area_nas'])
            full_open_dist = open_dist.merge(asymptotic_open_dist, on='term').merge(asymptotic_open_norm, on='term')
            print(full_open_dist)


            shut_dist = pd.DataFrame(shut_dist, columns=['term', 'w', 'rate', 'sh_tau_id', 'sh_area_id']).drop(['w', 'rate'], axis=1)
            asymptotic_shut_dist = pd.DataFrame(asymptotic_shut_dist, columns=['term', 'sh_tau_as', 'sh_area_as', 'rate']).drop(['rate'], axis=1)
            asymptotic_shut_norm = pd.DataFrame(asymptotic_shut_norm, columns=['term', 'sh_area_nas'])
            full_shut_dist = shut_dist.merge(asymptotic_shut_dist, on='term').merge(asymptotic_shut_norm, on='term')
            print(full_shut_dist)

            #'''
            #event_times_log = open(project + '_' + sc_type + '_' + file_name.strip('.abf') + '_event_times.txt', 'w')
            #n = event_times_log.write(event_times_data)
            #event_times_log.close()


            # THIS WORKS
            t, ipdf, epdf, apdf = scpl.shut_time_pdf(mec, sc_tres)
            plt.semilogx(t, ipdf, 'r--', t, epdf, 'b-', t, apdf, 'g-')
            plt.ylabel('fshut(t)')
            plt.xlabel('Shut time, ms')
            plt.title(sc_type + ' ' + file_name + ' ' + str(sc_tres*1000000) + ' ' + str(sc_tcrit*1000) )#+ '\n' + shuts_format)
            print('RED- ideal distribution\nGREEN- HJC distribution (corrected for missed events)')
            plt.savefig(project + '_' + sc_type + '_' + file_name.strip('.abf') + '_shut_plot.png')
            plt.close()
            plt.show()

            #THIS NOT
            t, ipdf, epdf, apdf = scpl.open_time_pdf(mec, sc_tres)
            plt.semilogx(t, ipdf, 'r--', t, epdf, 'b-', t, apdf, 'g-')
            plt.ylabel('fopen(t)')
            plt.xlabel('Open time, ms')
            plt.title(sc_type + ' ' + file_name + ' ' +  str(sc_tres * 1000000) + ' ' + str(
                sc_tcrit * 1000))  # + '\n' + shuts_format)
            print('RED- ideal distribution\nGREEN- HJC distribution (corrected for missed events)')
            plt.savefig(project + '_' + sc_type + '_' + file_name.strip('.abf') + '_open_plot.png')
            plt.close()
            plt.show()
            #'''
    # results saving below:

    if sc_model == 'CO':

        alpha = mec.Rates[0].rateconstants[0]
        beta = mec.Rates[1].rateconstants[0]

        refer_result = {'project': project, 'type': sc_type, 'file': file_name, 'model': sc_model,
                        'alpha': alpha, 'beta': beta,
                        }

    elif sc_model == 'CFO':

        alpha = mec.Rates[0].rateconstants[0]
        beta = mec.Rates[1].rateconstants[0]
        gamma = mec.Rates[2].rateconstants[0]
        delta = mec.Rates[3].rateconstants[0]

        refer_result = {'project': project, 'type': sc_type, 'file': file_name, 'model': sc_model,
                        'alpha': alpha, 'beta': beta,
                        'gamma': gamma, 'delta': delta,
                        't1_mod': t1_mod, 'p1_mod': p1_mod, 't2_mod': t2_mod, 'p2_mod': p2_mod,
                        }

    elif sc_model == 'CFOO':

        beta = mec.Rates[0].rateconstants[0]
        betap = mec.Rates[1].rateconstants[0]
        alpha = mec.Rates[2].rateconstants[0]
        alphap = mec.Rates[3].rateconstants[0]
        delta = mec.Rates[4].rateconstants[0]
        gamma = mec.Rates[5].rateconstants[0]

        refer_result = {'project': project, 'type': sc_type, 'file': file_name, 'model': sc_model,
                        'beta': beta, 'betap': betap, 'alpha': alpha, 'alphap': alphap, 'gamma': gamma}

    elif sc_model == 'CFOOD':

        beta = mec.Rates[0].rateconstants[0]
        betap = mec.Rates[1].rateconstants[0]
        alpha = mec.Rates[2].rateconstants[0]
        alphap = mec.Rates[3].rateconstants[0]
        delta = mec.Rates[4].rateconstants[0]
        gamma = mec.Rates[5].rateconstants[0]
        d = mec.Rates[6].rateconstants[0]
        r = mec.Rates[7].rateconstants[0]

        refer_result = {'project': project, 'type': sc_type, 'file': file_name, 'model': sc_model,
                        'beta': beta, 'betap': betap, 'alpha': alpha, 'alphap': alphap,
                        'gamma': gamma, 'delta': delta, 'd': d, 'r': r}

    elif sc_model == 'CFOD':

        beta = mec.Rates[0].rateconstants[0]
        alpha = mec.Rates[1].rateconstants[0]
        delta = mec.Rates[2].rateconstants[0]
        gamma = mec.Rates[3].rateconstants[0]
        d = mec.Rates[4].rateconstants[0]
        r = mec.Rates[5].rateconstants[0]

        refer_result = {'project': project, 'type': sc_type, 'file': file_name, 'model': sc_model,
                        'beta': beta, 'alpha': alpha,
                        'gamma': gamma, 'delta': delta, 'd': d, 'r': r}

    elif sc_model == 'CFODD':

        beta = mec.Rates[0].rateconstants[0]
        alpha = mec.Rates[1].rateconstants[0]
        delta = mec.Rates[2].rateconstants[0]
        gamma = mec.Rates[3].rateconstants[0]
        d = mec.Rates[4].rateconstants[0]
        dp = mec.Rates[5].rateconstants[0]
        r = mec.Rates[6].rateconstants[0]
        rp = mec.Rates[7].rateconstants[0]

        refer_result = {'project': project, 'type': sc_type, 'file': file_name, 'model': sc_model,
                        'delta': delta/1000, 'gamma': gamma/1000, 'beta': beta/1000,  'alpha': alpha/1000,
                        'd': d/1000, 'r': r/1000, 'dp': dp/1000, 'rp': rp/1000}

        model_result = {'project': project, 'type': sc_type, 'file': file_name, 'model': sc_model,
                        'delta': delta/1000, 'gamma': gamma/1000,
                        'beta': beta/1000, 'alpha': alpha/1000,
                        'd': d/1000, 'r': r/1000, 'dp': dp/1000, 'rp': rp/1000,
                        }
        full_open_long = full_open_dist.melt(id_vars='term', var_name='variable', value_name='value')
        full_open_long['var_term'] = full_open_long['variable'] + '_' + full_open_long['term'].astype(str)
        full_open_wide = full_open_long.pivot_table(index=None, columns='var_term', values='value')
        #print(df_wide)
        #full_open_dist_long = pd.DataFrame([df_wide.iloc[0]], columns=df_wide.columns)
        #print(full_open_dist_long)
        full_shut_long = full_shut_dist.melt(id_vars='term', var_name='variable', value_name='value')
        full_shut_long['var_term'] = full_shut_long['variable'] + '_' + full_shut_long['term'].astype(str)
        full_shut_wide = full_shut_long.pivot_table(index=None, columns='var_term', values='value')

        model_result = {**model_result, **full_open_wide.iloc[0].to_dict(), **full_shut_wide.iloc[0].to_dict()}


    elif sc_model == 'CFOODD':

        beta = mec.Rates[0].rateconstants[0]
        betap = mec.Rates[1].rateconstants[0]
        alpha = mec.Rates[2].rateconstants[0]
        alphap = mec.Rates[3].rateconstants[0]
        delta = mec.Rates[4].rateconstants[0]
        gamma = mec.Rates[5].rateconstants[0]
        d = mec.Rates[6].rateconstants[0]
        dp = mec.Rates[7].rateconstants[0]
        r = mec.Rates[8].rateconstants[0]
        rp = mec.Rates[9].rateconstants[0]

        refer_result = {'project': project, 'type': sc_type, 'file': file_name, 'model': sc_model,
                        'delta': delta/1000, 'gamma': gamma/1000,
                        'beta': beta/1000, 'alpha': alpha/1000, 'betap': betap/1000, 'alphap': alphap/1000,
                        'd': d/1000, 'r': r/1000, 'dp': dp/1000, 'rp': rp/1000}




    results.append(refer_result)
    model_results.append(model_result)

pd.DataFrame(model_results).to_csv('test_model_results.csv')

# REFER legacy
#results = pd.DataFrame(results)
#print(results)
#results.to_csv('hjcfit_rates_' + project + '.csv')


