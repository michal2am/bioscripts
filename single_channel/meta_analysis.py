import numpy as np
import pandas as pd
import plotly.express as px

meta = pd.read_csv('moje_meta_raw.csv', header=[0, 1])
models = pd.read_csv('moje_meta_models_raw.csv', header=[0, 1])
models.drop_duplicates(inplace=True, ignore_index=True)

merged = meta.merge(models, left_on=[('meta', 'file')], right_on=[('meta', 'file')], how='left')

merged.to_csv('moje_meta_merged_raw.csv')
merged.to_csv('moje_meta_merged_raw.csv')


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

# event time based REFER


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


merged = merged[merged[('meta', 'min_res')].notna()]                                                                    # only first row with full cell data

merged[('shuts', 't12_mean')] = merged.apply(lambda cell: mean_shut(cell), axis=1)                                      # calculate parameters


merged[('REFER_times', 'forward')] = merged.apply(lambda cell: forward(cell), axis=1)
print(merged[['meta','openings', 'shuts', 'REFER_times']])

merged[('REFER_times', 'equi')] = merged.apply(lambda cell: equi(cell), axis=1)
meta_REFERtimes = merged[['meta', 'REFER_times']].droplevel(0, axis=1)                                                  # flaten and select
meta_REFERtimes_grouped =meta_REFERtimes.groupby(by=['type', 'residue','residue_mut'])[['forward', 'equi']].mean()
# print(meta_REFERtimes_grouped)
meta_REFERtimes_grouped.reset_index(inplace=True)                                                                       # no multiindex for plotly
print(meta_REFERtimes_grouped)

for control, mutant in zip(['WT(F14/F31)', 'WT(F14/F31)', 'WT(F200)','WT(F45)', 'WT(F64)', 'WT(H55)', 'WT(P277)'],
                           ['F14', 'F31', 'F200', 'F45', 'F64', 'H55', 'P277']):

    cont_mut = meta_REFERtimes_grouped[(meta_REFERtimes_grouped['residue'] == mutant) | (meta_REFERtimes_grouped['type'] == control)]
    fig = px.scatter(cont_mut, x='equi', y='forward', title=mutant, hover_name='residue_mut')
    fig.add_trace(
        px.scatter(cont_mut, x='equi', y='forward',
                   trendline='ols',
                   color_discrete_sequence=px.colors.qualitative.Dark24,
                   ).data[1])

    fig.show()