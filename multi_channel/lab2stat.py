# parses channel lab results from csv file

import pandas as pd
from scipy.stats import ttest_ind

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def parse_model_rates():
    rates = pd.read_csv('simple_model_all_toCsv.csv')
    means = rates.groupby(['type']).mean()
    print(rates.round(2))
    print(means.round(2))
    means.to_csv('rates_means.csv')

    statistics = pd.DataFrame()

    for rec_type in ['L296V', 'G258V']:
        for rate in ['delta', 'gamma', 'beta', 'alpha', 'd', 'r']:

            wt = rates[rates.loc[:, 'type'] == 'WT'].loc[:, rate]
            mut = rates[rates.loc[:, 'type'] == rec_type].loc[:, rate]

            stat = ttest_ind(wt, mut)

            statistics = statistics.append(
                {'type': rec_type, 'rate': rate, 'p_value': stat.pvalue.round(6)}, ignore_index=True)

            if stat.pvalue < 0.05:
                print('Significant!')
                print('{}, {}, {}\n'.format(rate, rec_type, stat.pvalue.round(6)))

    for rec_type in ['L300V']:
        for rate in ['delta', 'gamma', 'beta', 'alpha', 'd', 'r', 'd', 'r']:

            wt = rates[rates.loc[:, 'type'] == 'WT_2D'].loc[:, rate]
            mut = rates[rates.loc[:, 'type'] == rec_type].loc[:, rate]

            stat = ttest_ind(wt, mut)

            statistics = statistics.append(
                {'type': rec_type, 'rate': rate, 'p_value': stat.pvalue.round(6)}, ignore_index=True)

            if stat.pvalue < 0.05:
                print('Significant!')
                print('{}, {}, {}\n'.format(rate, rec_type, stat.pvalue.round(6)))


    statistics.to_csv('rates_pvalues.csv')
    long_rates = rates.melt(id_vars=['type', 'cell'])
    return long_rates


def plot_model_rates(long_rates):

    sns.set_style()
    sns.set_context('talk')

    data = long_rates
    selected_types = ['WT', 'L296V', 'G258V']
    print(data[data.loc[:, 'type'] == 'G258V'])
    g = sns.catplot(
        data=data, kind="point",
        x='variable', y='value', hue='type',
        join=False, estimator=np.mean, ci=95,
        #order=selected_types,
        hue_order=selected_types,
        height=4, aspect=1.75,
        palette=sns.color_palette('deep'),
    )

    g.map(sns.swarmplot, 'variable', 'value', 'type',
          #order=selected_types,
          hue_order=selected_types,
          palette=sns.color_palette('deep'),
          alpha=0.5, marker='h')

    # g.axes[0, 0].axes.set_yticks(ticks=[0.5, 1, 2, 3])
    g.axes[0, 0].yaxis.set_major_locator(plt.MaxNLocator(5))

    g.despine(trim=True)
    g.set_axis_labels('rate', 'value [1/ms]')

    plt.savefig('rates_plot.png', dpi=600)
    plt.show()


parse_model_rates()
plot_model_rates(parse_model_rates())