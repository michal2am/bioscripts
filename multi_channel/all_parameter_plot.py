# plots all key parameters of macro traces

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

import seaborn as sns


selected_types = ['WT', 'L300V', 'L296V', 'G258V', 'WT_LC', 'G254V']
selected_who = ['II', 'AB', 'AD', 'AB']
selected_pulse = [500]

data = pd.read_csv('smarts_update_statistics_values_mm.csv', sep=';')
print(data)

data = data[data.loc[:, 'kto'].isin(selected_who) & data.loc[:, 'type'].isin(selected_types)]
data_long = pd.melt(data, id_vars=['kto', 'type', 'configuration', 'file', 'pulse_length', 'sweep'])

sns.set_style()
sns.set_context('talk')


def summary():
    summary = data.groupby(['type']).mean().loc[:,['rt_1090', 'FR10', 'FR300', 'FR500', 'des_Af', 'des_tf', 'des_As', 'des_ts', 'des_AC', 'dea_tm']].round(2)
    print(summary)

    statistics = pd.DataFrame()

    for parameter in ['rt_1090', 'FR10', 'FR300', 'FR500', 'des_Af', 'des_tf', 'des_As', 'des_ts', 'des_AC', 'dea_tm']:
        for rec_type in zip(['L300V', 'L296V', 'G258V', 'G254V'], ['WT', 'WT', 'WT', 'WT_LC']):

            wt = data[data.loc[:, 'type'] == rec_type[1]].loc[:, parameter]
            mut = data[data.loc[:, 'type'] == rec_type[0]].loc[:, parameter]

            stat = ttest_ind(wt, mut, nan_policy='omit')
            # print(stat)
            # print(wt, mut)
            if stat.pvalue < 0.05:
                print('Significant!')
                print('{}, {}, {}\n'.format(parameter, rec_type, stat.pvalue.round(6)))
                statistics = statistics.append({'type': rec_type[0], 'parameter': parameter, 'p_value': stat.pvalue.round(3)}, ignore_index=True)
            else:
                pass
                # print('Not significant!')
                # print('{}, {}, {}\n'.format(parameter, rec_type, stat.pvalue.round(6)))
    print(statistics)

    summary.to_csv('smarts_pandas_summary.csv')
    statistics.to_csv('smarts_pandas_statistics.csv')

def des_joint_plot():

    for des in [['des_Af', 'des_tf', 'des_f'], ['des_As', 'des_ts', 'des_s']]:

        g = sns.displot(data=data, kind='kde',
                        x=des[0], y=des[1], hue='type',
                        levels=1,
                        hue_order=selected_types,
                        height=4, aspect=1.75,
                        palette=sns.color_palette('deep', n_colors=len(selected_types)),
                        )
        g.map(sns.scatterplot, des[0], des[1], 'type',
              hue_order=selected_types,
              palette=sns.color_palette('deep', n_colors=len(selected_types)),
              )

        # g.axes[0, 0].axes.set_yticks(ticks=[0.5, 1, 2, 3])
        g.axes[0, 0].yaxis.set_major_locator(plt.MaxNLocator(5))
        g.axes[0, 0].xaxis.set_major_locator(plt.MaxNLocator(5))

        g.despine(trim=True)
        g.set_axis_labels(des[0], des[1])

        plt.savefig(des[2] + '.png', dpi=600)
        plt.show()


def interval_plot(param):

    g = sns.catplot(
        data=data, kind="point",
        x="type", y=param, hue='type',
        join=False, estimator=np.mean, ci=95,
        order=selected_types,
        hue_order=selected_types,
        height=4, aspect=1.75,
        palette=sns.color_palette('deep'),
    )

    g.map(sns.swarmplot, "type", param, 'type',
          order=selected_types,
          hue_order=selected_types,
          palette=sns.color_palette('deep'),
          alpha=0.5, marker='h')

    # g.axes[0, 0].axes.set_yticks(ticks=[0.5, 1, 2, 3])
    g.axes[0, 0].yaxis.set_major_locator(plt.MaxNLocator(5))

    g.despine(trim=True)
    g.set_axis_labels("", param)

    plt.savefig(param + '.png', dpi=600)
    plt.show()


def fr_plot():

    g = sns.catplot(
        data=data_long[data_long['variable'].isin(['FR10', 'FR300', 'FR500'])], kind='point',
        x="variable", y="value", hue='type',
        join=True, dodge=False, estimator=np.mean, ci=95,
        order=['FR10', 'FR300', 'FR500'],
        hue_order=selected_types,
        height=4, aspect=1.75,
        #palette=sns.color_palette('deep', n_colors=len(selected_types)),
    )

    g.map(sns.swarmplot, 'variable', 'value', 'type',
          order=['FR10', 'FR300', 'FR500'],
          hue_order=selected_types,
          #palette=sns.color_palette('deep', n_colors=len(selected_types)),
          alpha=0.5, marker='h')

    g.axes[0, 0].axes.set_yticks(ticks=[0, .25, .5, .75, 1])
    #g.axes[0, 0].yaxis.set_major_locator(plt.MaxNLocator(5))

    g.despine(trim=True)
    g.set_axis_labels("", "")
    g.legend.set_title("receptor type")

    plt.savefig('fr_all.png', dpi=600)
    plt.show()


def des_a_plot():

    g = sns.catplot(
        data=data_long[data_long['variable'].isin(['des_Af', 'des_As', 'des_AC'])], kind='point',
        x="variable", y="value", hue='type',
        join=True, dodge=False, estimator=np.mean, ci=95,
        order=['des_Af', 'des_As', 'des_AC'],
        hue_order=selected_types,
        height=4, aspect=1.75,
        palette=sns.color_palette('deep'),
    )

    g.map(sns.swarmplot, 'variable', 'value', 'type',
          order=['des_Af', 'des_As', 'des_AC'],
          hue_order=selected_types,
          palette=sns.color_palette('deep'),
          alpha=0.5, marker='h')

    g.axes[0, 0].axes.set_yticks(ticks=[0, .25, .5, .75, 1])
    #g.axes[0, 0].yaxis.set_major_locator(plt.MaxNLocator(5))

    g.despine(trim=True)
    g.set_axis_labels("", "")
    g.legend.set_title("receptor type")

    plt.savefig('des_a_all.png', dpi=600)
    plt.show()


summary()

plots = False

if plots:

    interval_plot('rt_1090')
    interval_plot('dea_tm')
    interval_plot('des_tf')
    interval_plot('des_ts')

    fr_plot()

    des_joint_plot()
    des_a_plot()
