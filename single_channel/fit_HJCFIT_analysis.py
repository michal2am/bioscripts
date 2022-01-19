import pandas as pd
import plotly.express as px
import numpy as np
from scikit_posthocs import outliers_iqr


full_data = pd.read_csv('moje_tcrits_raw_iter.csv', sep=',')
full_data['tcrit'] = full_data['tcrit'].astype('str')
#print(full_data)

full_data_clear = pd.DataFrame(columns=full_data.columns)

for mut in full_data.type.unique():
    for tcrit in full_data.tcrit.unique():
        for param in ['alpha', 'beta', 'delta', 'gamma', 't1_mod', 't2_mod']:

            selected = full_data[(full_data['type'] == mut) & (full_data['tcrit'] == tcrit)].copy()

            outliers_val = outliers_iqr(selected[param], ret='outliers')
            print(mut, tcrit, param, outliers_val)
            full_data[param].replace(outliers_val, np.NaN, inplace=True)
            #selected[param].replace(outliers_val, np.NaN, inplace=True)
            #print(selected[param])
            #full_data_clear = full_data_clear.append(selected, ignore_index=True)

#full_data_clear.to_csv('test.csv')

#print(full_data_clear['alpha'])
#print(full_data_clear['beta'])

for rate in ['alpha', 'beta', 'gamma', 'delta', 't1_mod', 't2_mod']:
    pass
    # fig = px.bar(full_data, x='file', y=rate, barmode='group', color='tcrit',hover_name='type', title='xxx')
    #fig.show()
    fig2 = px.box(full_data, x='type', y=rate, color='tcrit', hover_name='type', title='xxx')
    fig2.show()

