import numpy as np
import pandas as pd
import plotly.express as px



def rates_vs_res():

    model_vs_res = merged[['rates', 'meta']].droplevel(0, axis=1)
    model_vs_res = model_vs_res[model_vs_res['min_res'].notna()]                                                        # only first row with full cell data
    model_vs_res['min_res'] = pd.to_numeric(model_vs_res['min_res'])                                                    # res as numeric

    print(model_vs_res)

    for rate in ['b2', 'a2', 'gamma', 'delta']:
        fig = px.scatter(model_vs_res, x='min_res', y=rate, color='residue_mut', symbol='type')
        fig.show()
        # fig = px.scatter(model_vs_res, x='fc', y=rate, color='residue_mut', symbol='type')
        # fig.show()


def REFER_time(merged):
    """event time based REFER"""

    def mean_shut(cell):
        percent_p1_p2 = cell[('shuts', '%P1')] + cell[('shuts', '%P2')]
        norm_P1 = cell[('shuts', '%P1')] / percent_p1_p2
        norm_P2 = cell[('shuts', '%P2')] / percent_p1_p2
        mean_t1_t2 = cell[('shuts', 't1')] * norm_P1 + cell[('shuts', 't2')] * norm_P2
        return mean_t1_t2

    def forward(cell):
        return np.log(1/cell[('shuts', 't12_mean')])

    def equi(cell):
        return np.log(cell[('openings', 't_mean')]/cell[('shuts', 't12_mean')])

    merged = merged[merged[('meta', 'min_res')].notna()]                                                                # only first row with full cell data

    merged[('shuts', 't12_mean')] = merged.apply(lambda cell: mean_shut(cell), axis=1)                                  # calculate parameters
    merged[('REFER_times', 'forward')] = merged.apply(lambda cell: forward(cell), axis=1)
    merged[('REFER_times', 'equi')] = merged.apply(lambda cell: equi(cell), axis=1)

    meta_REFERtimes = merged[['meta', 'REFER_times']].droplevel(0, axis=1)                                              # flaten and select
    meta_REFERtimes_grouped =meta_REFERtimes.groupby(by=['type', 'residue','residue_mut'])[['forward', 'equi']].mean()
    # print(meta_REFERtimes_grouped)

    meta_REFERtimes_grouped.reset_index(inplace=True)                                                                   # no multiindex for plotly
    print(meta_REFERtimes_grouped)

    for control, mutant in zip(['WT(F14/F31)', 'WT(F14/F31)', 'WT(F200)','WT(F45)', 'WT(F64)', 'WT(H55)', 'WT(P277)'],
                               ['F14', 'F31', 'F200', 'F45', 'F64', 'H55', 'P277']):

        cont_mut = meta_REFERtimes_grouped[(meta_REFERtimes_grouped['residue'] == mutant) | (meta_REFERtimes_grouped['type'] == control)]
        fig = px.scatter(cont_mut, x='equi', y='forward', title="{} REFER by event time".format(mutant),
                         hover_name='residue_mut', color='residue_mut',
                         template='presentation', width=600, height=600,)
        fig.add_trace(
            px.scatter(cont_mut, x='equi', y='forward',
                       trendline='ols',
                       color_discrete_sequence=px.colors.qualitative.Dark24,
                       ).data[1])
        fig.show()

        with open(mutant + '_auerbach_times.html', 'w') as f:
            f.write(fig.to_html())


def prepare_hjcfit_config(mutant, control, file_name, model, tcrit):
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

    selected = merged[(merged[('meta', 'residue')] == mutant) | (merged[('meta', 'type')] == control)]

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
    # print(selected)

    selected.to_csv(file_name)


# Bambi's meta file with all the data
meta = pd.read_csv('moje_meta_raw.csv', header=[0, 1])
# my file with Bambi's modeling results, just a copy with rows corresponding to the cell file name and rates
models = pd.read_csv('moje_meta_models_raw.csv', header=[0, 1])
models.drop_duplicates(inplace=True, ignore_index=True)                                                                 # necessary for duplicated WT cells

merged = meta.merge(models, left_on=[('meta', 'file')], right_on=[('meta', 'file')], how='left')                        # basically all Bambi's spreadsheets in one dataframe
merged.to_csv('moje_meta_merged_raw.csv')


# PLAYGROUND BELOW

REFER_time(merged)

prepare_hjcfit_config('F200', 'WT(F200)', 'hjcfit_config_f200_MetaBambiCO.csv', 'CO', 'tcrit_CFO')
prepare_hjcfit_config('F64', 'WT(F64)', 'hjcfit_config_f64_MetaBambiCO.csv', 'CO', 'tcrit_CFO')     # f64g?
prepare_hjcfit_config('P277', 'WT(P277)', 'hjcfit_config_p277_MetaBambiCO.csv', 'CO', 'tcrit_CFO')
prepare_hjcfit_config('F45', 'WT(F45)', 'hjcfit_config_f45_MetaBambiCO.csv', 'CO', 'tcrit_CFO')
prepare_hjcfit_config('F14', 'WT(F14/F31)', 'hjcfit_config_f14_MetaBambiCO.csv', 'CO', 'tcrit_CFO')
prepare_hjcfit_config('F31', 'WT(F14/F31)', 'hjcfit_config_f31_MetaBambiCO.csv', 'CO', 'tcrit_CFO')
prepare_hjcfit_config('H55', 'WT(H55)', 'hjcfit_config_h55_MetaBambiCO.csv', 'CO', 'tcrit_CFO')
