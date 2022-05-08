import pandas as pd


def prepare_hjcfit_config(meta, mutant, control, file_name, model, tcrit):
    '''
    parses meta-file with Bambi analysis for input into python-hjcfit
    TODO: supports only CO-model
    TODO: no support for exp t_x insertion into config
    :param mutant: e.g. 'F200'
    :param control: e.g 'WT(F200)'
    :param file_name: e.g 'hjcfit_config_f200_MetaBambiCO.csv'
    :param model: e.g 'CO'
    :param tcrit: tcrit_CFO (first KT approx) or final_tcrit (second KT approx)
    :return: csv file ready for hjcfit
    '''

    selected = meta[(meta[('meta', 'residue')] == mutant) | (meta[('meta', 'type')] == control)]

    selected = selected.loc[:, ['meta', 't_crit']].droplevel(0, axis=1)
    selected = selected.loc[:, ['residue_mut', 'file', 'cluster_name', 'min_res', tcrit]]

    # careful here! missing data may be copied from the wrong cell, works fine only for copying info for subsequent scns from same cell
    selected.fillna(method='ffill', inplace=True)
    print('Cells with t critical 0, will be removed ... ')
    print(selected[selected[tcrit] == 0])
    selected.drop(selected[selected[tcrit] == 0].index, inplace=True)

    selected['beta'] = 5000
    selected['alpha'] = 5000
    selected['model'] = model

    selected.rename(columns={'residue_mut': 'type', 'cluster_name': 'file_scn', 'min_res': 'tres', tcrit: 'tcrit'}, inplace=True)

    selected.replace({'brak': 'WT'}, inplace=True)
    print(selected)

    selected.to_csv(file_name)


# Bambi's meta file with all the data
meta = pd.read_csv('moje_meta_raw.csv', header=[0, 1])
meta.to_csv('moje_meta_merged_raw.csv')

#prepare_hjcfit_config('F200', 'WT(F200)', 'hjcfit_config_f200_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('F64', 'WT(F64)', 'hjcfit_config_f64_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('P277', 'WT(P277)', 'hjcfit_config_p277_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('F45', 'WT(F45)', 'hjcfit_config_f45_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('F14', 'WT(F14/F31)', 'hjcfit_config_f14_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('F31', 'WT(F14/F31)', 'hjcfit_config_f31_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('H55', 'WT(H55)', 'hjcfit_config_h55_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('L296', 'WT(F14/F31)', 'hjcfit_config_l296_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('E270', 'WT(F14/F31)', 'hjcfit_config_e270_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('H267', 'WT(F14/F31)', 'hjcfit_config_e270_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('P273', 'WT(F45)', 'hjcfit_config_p273_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('E153', 'WT(E153)', 'hjcfit_config_e153_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('F45', 'WT(E153)', 'hjcfit_config_f45wte153_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('F45', 'WT(F14/F31)', 'hjcfit_config_f45wtBambi_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('F200', 'WT(F14/F31)', 'hjcfit_config_f200wtBambi_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('F200', 'WT(E153)', 'hjcfit_config_f200wte153_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
#prepare_hjcfit_config('E270', 'WT(F14/F31)', 'hjcfit_config_e270_MetaBambiCOfina.csv', 'CO', 'final_tcrit')
prepare_hjcfit_config(meta, 'V53', 'WT(V53)', 'hjcfit_config_v53_MetaBambiCOfina.csv', 'CO', 'final_tcrit')