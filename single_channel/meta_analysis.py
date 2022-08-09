import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scikit_posthocs import outliers_iqr


sns.set_style()
sns.set_context("talk")
sns.set_palette('muted')


def outliers_iqrsckit(data, rates_list):
    for mut in data.meta_residue_mut.unique():
        print('Looking for outliers in {}'.format(mut))
        for rate in rates_list:
            print('Checking rate {}'.format(rate))
            rates = list(
                data.loc[data['meta_residue_mut'] == mut][rate])  # list() just for index reset
            cells = list(data.loc[data['meta_residue_mut'] == mut]['meta_file'])
            outliers_val = outliers_iqr(rates, ret='outliers')
            outliers_ind = outliers_iqr(rates, ret='outliers_indices')

            if outliers_val.size > 0:
                outliers_cel = [cells[ind] for ind in outliers_ind]
                print('{} value outliers in {} type; cells: {}, values: {} \n within: {}'.
                      format(rate, mut, outliers_cel, outliers_val, rates))
                data.loc[data['meta_file'].isin(outliers_cel), rate] = np.NAN

meta = pd.read_csv('moje_meta_raw_WTMonikaOut.csv', header=[0, 1])
#print(meta[['meta', 'openings', 'shuts']])

meta.columns = ['_'.join(col) for col in meta.columns.values]
meta.dropna(subset=['meta_cells_no'], inplace=True)

#meta = meta[~meta['meta_residue'].isin(['F14'])]
#meta = meta[~meta['meta_residue_mut'].isin(['F64C', 'F64L', 'V53E'])]

#print(meta)

rates = (meta[['meta_file', 'meta_residue', 'meta_residue_mut', 'shuts_t1', 'shuts_t2', 'shuts_t3', 'shuts_t4', 'shuts_%P1', 'shuts_%P2', 'shuts_%P3', 'shuts_%P4',
               'openings_t1', 'openings_t2',  'openings_%P1', 'openings_%P2']])
#cumulative_system = datasets.groupby(['System', 'State', 'Lipid'], as_index=False)['Event Threshold'].sum()

print(rates)
rates_list = ['shuts_t1', 'shuts_t2', 'shuts_t3', 'shuts_t4', 'shuts_%P1', 'shuts_%P2', 'shuts_%P3', 'shuts_%P4',
               'openings_t1', 'openings_t2',  'openings_%P1', 'openings_%P2']
outliers_iqrsckit(rates, rates_list)
print(rates)

rates = rates.groupby(['meta_file','meta_residue', 'meta_residue_mut'], as_index=False).mean()

correlations = rates.corr()
mask = np.zeros_like(correlations)
mask[np.triu_indices_from(mask)] = True
correlations[(correlations >= -.25) & (correlations <= .25)] = np.nan

fig, ax = plt.subplots(figsize=(15,15))

#sns.heatmap(correlations, annot=True, mask=mask)
sns.heatmap(correlations, annot=True)

plt.savefig('correlations' + '_sns.png', dpi=300)

fig, ax = plt.subplots(figsize=(15,15))

g = sns.pairplot(rates, diag_kind='kde', hue='meta_residue')
plt.savefig('raw_correlations' + '_sns.png', dpi=300)
#plt.show()
