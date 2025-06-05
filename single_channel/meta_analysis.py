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

meta = pd.read_csv('REFER_zbiorczy_2024_raw.csv', header=[0, 1])

meta.columns = ['_'.join(col) for col in meta.columns.values]
meta.dropna(subset=['meta_cells_no'], inplace=True)




rates = (meta[['meta_file', 'meta_type','meta_residue', 'meta_residue_mut', 'meta_min_res', 'shuts_t1', 'shuts_t2', 'shuts_t3', 'shuts_t4', 'shuts_%P1', 'shuts_%P2', 'shuts_%P3', 'shuts_%P4',
               'openings_t1', 'openings_t2',  'openings_%P1', 'openings_%P2']])

rates_list = ['meta_min_res', 'shuts_t1', 'shuts_t2', 'shuts_t3', 'shuts_t4', 'shuts_%P1', 'shuts_%P2', 'shuts_%P3', 'shuts_%P4',
               'openings_t1', 'openings_t2',  'openings_%P1', 'openings_%P2']
#outliers_iqrsckit(rates, rates_list)
rates.to_csv('test')
print(rates)

#rates = rates.groupby(['meta_file','meta_residue', 'meta_residue_mut'], as_index=False).mean()

#print(rates)


rates_WT = rates[rates['meta_residue'] == 'brak']
print(rates_WT)
#rates_WT.to_csv('wt_manual.csv')
#sns.scatterplot(x=rates_WT['openings_t1'], y=rates_WT['openings_%P2'])
#plt.show()

def plot(data, feature, yticks=False):
    sns.set_style()
    sns.set_context("paper")

    g = sns.catplot(
        data=data, kind="point",
        x="meta_type", y=feature,  join=False, estimator=np.mean, errorbar=('ci', 98),
        height=3, aspect=3,
    )

    g.map(sns.swarmplot, "meta_type", feature, alpha=.2,)

    if yticks:
        g.axes[0, 0].axes.set_yticks(ticks=yticks)

    g.despine(trim=True)
    g.set_axis_labels("", feature)
    plt.savefig(feature + '.png', dpi=300)
    plt.show()

#plot(rates_WT,'openings_%P1', [0, 0.25, 0.5, 0.75, 1])
#plot(rates_WT,'openings_%P2', [0, 0.25, 0.5, 0.75, 1])
#plot(rates_WT,'openings_t1', [0, 1, 2, 3, 4, 5])
#plot(rates_WT,'openings_t2', [0, 1, 2, 3, 4, 5, 7, 8 ,9, 10])

corr = True
if corr:
    correlations = rates_WT[rates_list].corr()
    mask = np.zeros_like(correlations)
    mask[np.triu_indices_from(mask)] = True
    correlations[(correlations >= -.25) & (correlations <= .25)] = np.nan

    fig, ax = plt.subplots(figsize=(15,15))

    #sns.heatmap(correlations, annot=True, mask=mask)
    sns.heatmap(correlations, annot=True)

    plt.savefig('correlationsWT' + '_sns.png', dpi=300)

    fig, ax = plt.subplots(figsize=(15,15))

    g = sns.pairplot(rates_WT, diag_kind='kde', hue='meta_residue')
    plt.savefig('raw_correlationsWT' + '_sns.png', dpi=300)
    #plt.show()
