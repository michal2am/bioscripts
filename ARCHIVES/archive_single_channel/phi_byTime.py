import numpy as np
import pandas as pd
import plotly.express as px
from scikit_posthocs import outliers_iqr

def REFER_time(merged):
    """event time based REFER"""

    def mean_shut(cell):
        percent_p1_p2 = cell[('shuts', '%P1')] + cell[('shuts', '%P2')]
        norm_P1 = cell[('shuts', '%P1')] / percent_p1_p2
        norm_P2 = cell[('shuts', '%P2')] / percent_p1_p2
        mean_t1_t2 = cell[('shuts', 't1')] * norm_P1 + cell[('shuts', 't2')] * norm_P2
        return mean_t1_t2

    def alpha(cell):
        return 1/cell[('openings', 't_mean')]

    def beta(cell):
        return 1/cell[('shuts', 't12_mean')]

    def forward(cell):
        return np.log(cell['beta'])

    def equi(cell):
        return np.log(cell['beta']/cell['alpha'])

    def outliers_iqrsckit(flat_cells):

        for mut in flat_cells.residue_mut.unique():
            # TODO: outliers for WT
            if mut == 'brak': continue
            print('Looking for outliers in {}'.format(mut))
            for rate in ['alpha', 'beta']:
                print('Checking rate {}'.format(rate))
                rates = list(
                    flat_cells[flat_cells['residue_mut'] == mut][rate])  # list() just for index reset
                cells = list(flat_cells[flat_cells['residue_mut'] == mut]['file'])
                outliers_val = outliers_iqr(rates, ret='outliers')
                outliers_ind = outliers_iqr(rates, ret='outliers_indices')

                if outliers_val.size > 0:
                    outliers_cel = [cells[ind] for ind in outliers_ind]
                    print('{} value outliers in {} type; cells: {}, values: {} \n within: {}'.
                          format(rate, mut, outliers_cel, outliers_val, rates))
                    flat_cells.loc[flat_cells['file'].isin(outliers_cel), 'alpha'] = np.NAN

    # TODO: criterium for cell selection
    pd.set_option('display.max_rows', None)
    # merged = merged[merged[('meta', 'min_res')].notna()]                                                              # only first row with full cell data
    merged = merged[merged[('meta', 'clusters_no')].notna()]                                                            # only first row with full cell data

    merged[('shuts', 't12_mean')] = merged.apply(lambda cell: mean_shut(cell), axis=1)                                  # calculate parameters
    merged[('REFER_times', 'alpha')] = merged.apply(lambda cell: alpha(cell), axis=1)
    merged[('REFER_times', 'beta')] = merged.apply(lambda cell: beta(cell), axis=1)

    #print(merged[['meta', 'REFER_times']])

    meta_REFERtimes = merged[['meta', 'REFER_times']].droplevel(0, axis=1)                                              # flaten and select for outliers and plots
    outliers_iqrsckit(meta_REFERtimes)

    meta_REFERtimes['forward'] = meta_REFERtimes.apply(lambda cell: forward(cell), axis=1)
    meta_REFERtimes['equi'] = meta_REFERtimes.apply(lambda cell: equi(cell), axis=1)
    # print(meta_REFERtimes)

    cell_averages = True
    if cell_averages:
        meta_REFERtimes = meta_REFERtimes.groupby(by=['type', 'residue','residue_mut'])[['forward', 'equi', 'alpha', 'beta']].mean()

    meta_REFERtimes.reset_index(inplace=True)  # no multiindex for plotly
    print(meta_REFERtimes)

    controls = ['WT(F14/F31)', 'WT(F14/F31)', 'WT(F200)','WT(F45)', 'WT(F64)', 'WT(H55)', 'WT(P277)', 'WT(F14/F31)', 'WT(F14/F31)', 'WT(F45)', 'WT(F14/F31)', 'WT(F14/F31)', 'WT(E153)', 'WT(E153)']
    mutants = ['F14', 'F31', 'F200', 'F45', 'F64', 'H55', 'P277', 'L296', 'L300', 'P273', 'H267', 'E270', 'E153', 'V53']
    # controls = ['WT(F14/F31)', 'WT(F14/F31)', 'WT(F14/F31)']
    # mutants = ['H267', 'L296', 'L300',]

    for control, mutant in zip(controls, mutants):

        print(control, mutant)

        cont_mut = meta_REFERtimes[(meta_REFERtimes['residue'] == mutant) | (meta_REFERtimes['type'] == control)]

        plot = (px.scatter(cont_mut, x='equi', y='forward',
                           # title="{} REFER by HJCFIT CO, WT from {}".format(project, WT_project),
                           title="{}".format(mutant.upper()),
                           labels={'equilibrium_raw': 'log(equilibrium rate)',
                                   'forward_raw': 'log(forward rate)'},
                           color='type', template='presentation', width=400, height=400,
                           hover_name='type', hover_data=['alpha', 'beta'],
                           color_discrete_sequence=px.colors.qualitative.Dark24,
                           ))
        plot.add_trace(
            px.scatter(cont_mut, x='equi', y='forward',
                       trendline='ols',
                       color_discrete_sequence=px.colors.qualitative.Dark24,
                       ).data[1])

        plot.write_html(mutant + '_auerbach_times.html',)

# Bambi's meta file with all the data
meta = pd.read_csv('moje_meta_raw.csv', header=[0, 1])

REFER_time(meta)
