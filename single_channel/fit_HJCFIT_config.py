import pandas as pd
import numpy as np
import argparse
from scipy.optimize import fsolve


def prepare_hjcfit_config(meta, mutant, control, file_name, model):
    '''
    parses meta-file with Bambi analysis for input into python-hjcfit
    TODO: supports only CO-model
    TODO: no support for exp t_x insertion into config
    :param mutant: e.g. 'F200'
    :param control: e.g 'WT(F200)'
    :param file_name: e.g 'hjcfit_config_f200_MetaBambiCO.csv'
    :param model: e.g 'CO'
    :return: csv file ready for hjcfit
    '''

    selected = meta[(meta[('meta', 'residue')] == mutant) | (meta[('meta', 'type')] == control)]

    selected = selected.loc[:, ['meta', 'shuts', 't_crit']].droplevel(0, axis=1)
    selected = selected.loc[:, ['residue_mut', 'file', 'cluster_name', 'min_res', 'final_tcrit', 't2', 't3', '%P2', '%P3']]

    # careful here! missing data may be copied from the wrong cell, works fine only for copying info for subsequent scns from same cell
    selected.fillna(method='ffill', inplace=True)
    print('Cells with t critical 0, will be removed ... ')
    print(selected[selected['final_tcrit'] == 0])
    selected.drop(selected[selected['final_tcrit'] == 0].index, inplace=True)

    selected['model'] = model

    if model == 'CO':
        selected['beta'] = 5000
        selected['alpha'] = 5000

    if model == 'CFOODD':
        selected['beta'] = 15000
        selected['betap'] = 5000
        selected['alpha'] = 1000
        selected['alphap'] = 1000
        selected['delta'] = 5000
        selected['gamma'] = 5000
        selected['d'] = 1000
        selected['dp'] = 500
        selected['r'] = 100
        selected['rp'] = 50

    if model == 'CFODD':
        selected['beta'] = 15000
        selected['alpha'] = 1000
        selected['delta'] = 5000
        selected['gamma'] = 5000
        selected['d'] = 1000
        selected['dp'] = 500
        selected['r'] = 100
        selected['rp'] = 50

    selected.rename(columns={'residue_mut': 'type', 'cluster_name': 'file_scn', 'min_res': 'tres',
                             'final_tcrit': 'tcrit_KT'}, inplace=True)

    selected.replace({'brak': 'WT'}, inplace=True)

    selected['tcrit_CS'] = selected.apply(calculate_tcrit_CS, axis=1)
    selected['tcrit_J'] = selected.apply(calculate_tcrit_J, axis=1)
    selected['tcrit_MPCN'] = selected.apply(calculate_tcrit_MPCN, axis=1)
    selected['tcrit_inf'] = 1000000

    print(selected)
    selected.to_csv(file_name)


def calculate_tcrit_CS(config_row):
    tcrit_guess = (config_row['t2'] + config_row['t3'])/2
    func_tcrit_DC = lambda tcrit_DC : 1 - np.exp(-tcrit_DC/config_row['t3']) - np.exp(-tcrit_DC/config_row['t2'])
    tcrit_DC_solution = fsolve(func_tcrit_DC, tcrit_guess)
    return tcrit_DC_solution[0]

def calculate_tcrit_MPCN(config_row):
    tcrit_guess = (config_row['t2'] + config_row['t3'])/2
    func_tcrit_MPCN = lambda tcrit_MPCN : config_row['%P3']*(1 - np.exp(-tcrit_MPCN/config_row['t3'])) - config_row['%P2']*np.exp(-tcrit_MPCN/config_row['t2'])
    tcrit_MPCN_solution = fsolve(func_tcrit_MPCN, tcrit_guess)
    return tcrit_MPCN_solution[0]

def calculate_tcrit_J(config_row):
    tcrit_guess = (config_row['t2'] + config_row['t3'])/2
    func_tcrit_J = lambda tcrit_J : (config_row['%P2']/config_row['t2'])*np.exp(-tcrit_J/config_row['t2']) - (config_row['%P3']/config_row['t3'])*np.exp(-tcrit_J/config_row['t3'])
    tcrit_J_solution = fsolve(func_tcrit_J, tcrit_guess)
    return tcrit_J_solution[0]


parser = argparse.ArgumentParser()
parser.add_argument('--meta_file')

parser.add_argument('--mutant')
parser.add_argument('--control')
parser.add_argument('--file_name')
parser.add_argument('--model')

args = parser.parse_args()

# Bambi's meta file with all the data
meta = pd.read_csv(args.meta_file, header=[0, 1])
print(meta)

prepare_hjcfit_config(meta, args.mutant, args.control, args.file_name, args.model)

#prepare_hjcfit_config(meta, 'F45', 'WT(F45)', 'hjcfit_config_f45_2024_tcritm.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config(meta, 'F200', 'WT(F200)', 'hjcfit_config_f200_2024_4tcrits.csv', 'CO')
#prepare_hjcfit_config(meta, 'F31', 'WT(F14/F31)', 'hjcfit_config_f31_2024_tcritm.csv', 'CO', 'final_tcrit')


#prepare_hjcfit_config(meta, 'V53', 'WT(V53)', 'hjcfit_config_v53_fullModels.csv', 'CFOODD', 'inf_tcrit')

#prepare_hjcfit_config(meta, 'F14', 'WT(F14/F31)', 'hjcfit_config_f14_2022.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config(meta, 'F31', 'WT(F14/F31)', 'hjcfit_config_f31_2022.csv', 'CO', 'final_tcrit')

#prepare_hjcfit_config(meta, 'F200', 'WT(F200)', 'hjcfit_config_f200_2022.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config(meta, 'E153', 'WT(E153)', 'hjcfit_config_e153_2022.csv', 'CO', 'final_tcrit')

#prepare_hjcfit_config(meta, 'F64', 'WT(F64)', 'hjcfit_config_f64_2022.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config(meta, 'F45', 'WT(F45)', 'hjcfit_config_f45_2022.csv', 'CO', 'final_tcrit')

#prepare_hjcfit_config(meta, 'V53', 'WT(V53)', 'hjcfit_config_v53_2022.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config(meta, 'H55', 'WT(H55)', 'hjcfit_config_h55_2022a.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config(meta, 'P277', 'WT(P277)', 'hjcfit_config_p277_2022.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config(meta, 'P273', 'WT(F45)', 'hjcfit_config_p273_2022.csv', 'CO', 'final_tcrit')

#prepare_hjcfit_config(meta, 'E270', 'WT(E153)', 'hjcfit_config_e270_2022.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config(meta, 'H267', 'WT(E153)', 'hjcfit_config_h267_2022.csv', 'CO', 'final_tcrit')

#prepare_hjcfit_config(meta, 'L300', 'WT(E153)', 'hjcfit_config_l300_2022.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config(meta, 'L296', 'WT(E153)', 'hjcfit_config_l296_2022.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config(meta, 'G258', 'WT(G258)', 'hjcfit_config_g258_2022.csv', 'CO', 'final_tcrit')
