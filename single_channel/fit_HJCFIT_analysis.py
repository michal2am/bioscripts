import pandas as pd
import plotly.express as px
import numpy as np
from scikit_posthocs import outliers_iqr


full_data = pd.read_csv('moje_tcrits_raw.csv', sep=';')
full_data['tcrit'] = full_data['tcrit'].astype('str')
print(full_data)

full_data_clear = pd.DataFrame()

for mut in full_data.type.unique():
    for tcrit in full_data.tcrit.unique():
        for param in ['alpha', 'beta', 'delta', 'gamma']:

            selected = full_data[(full_data['type'] == mut) & (full_data['tcrit'] == tcrit)].copy()

            outliers_val = outliers_iqr(selected[param], ret='outliers')
            selected[param].replace(outliers_val, np.NaN, inplace=True)
            print(selected)
            full_data_clear.append(selected, ignore_index=True)

print(full_data_clear)
print(full_data_clear['beta'])

for rate in ['alpha', 'beta', 'gamma', 'delta', 't1_mod', 't2_mod']:
    fig = px.bar(full_data, x='file', y=rate, barmode='group', color='tcrit',hover_name='type', title='xxx')
    #fig.show()
    fig2 = px.box(full_data, x='type', y=rate, color='tcrit', hover_name='type', title='xxx')
    #fig2.show()

'''
def clear_outliers(data, params, category):

    for mut in data.type.unique():
        for tcrit in ['50', '75', '100']:
            for param in params:

                rates = list(data[(data['type'] == mut) & (data['tcrit'] == tcrit)][param])  # list() just for index reset
                #cells = list(self.table_results.loc[self.table_results['type'] == mut]['file'])
                outliers_val = outliers_iqr(rates, ret='outliers')
                outliers_ind = outliers_iqr(rates, ret='outliers_indices')

                if outliers_val.size > 0:
                    outliers_cel = [cells[ind] for ind in outliers_ind]
                    print('{} value outliers in {} type; cells: {}, values: {} \n within: {}'.
                          format(rate, mut, outliers_cel, outliers_val, rates))
                    self.table_results.loc[self.table_results['file'].isin(outliers_cel), rate] = np.NAN
'''