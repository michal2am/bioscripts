import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import f_oneway
from scipy.stats import kruskal
from scipy.stats import shapiro
from scipy.stats import levene
from scikit_posthocs import outliers_iqr
from scikit_posthocs import posthoc_dunn
from scikit_posthocs import posthoc_tukey
from scikit_posthocs import posthoc_dunnett
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels import robust

import argparse
import matplotlib.pyplot as plt

sns.set_style()
sns.set_context("talk")
sns.set_palette('muted')


def outliers_mz(data, rates_list):

    def modified_zscore(series):
        median = series.median()
        mad = np.median(np.abs(series - median))
        if mad == 0:
            return pd.Series(np.zeros(len(series)), index=series.index)
        return 0.6745 * (series - median) / mad

    for rate in rates_list:
        z_scores = data.groupby('meta_receptor')[rate].transform(modified_zscore)
        outlier_mask = np.abs(z_scores) > 3.5
        data[f'{rate}_outlier'] = outlier_mask

        print(f'\n=== Analysis for "{rate}" ===')
        for group, group_df in data.groupby('meta_receptor'):
            original_values = group_df[rate]
            outliers = original_values[group_df[f'{rate}_outlier']]
            print(f'\nGroup: {group}')
            print(f'  All values: {original_values.tolist()}')
            print(f'  Outliers:   {outliers.tolist()}')

        data.loc[outlier_mask, rate] = np.NaN



def outliers_iqrsckit(data, rates_list):
    for mut in data.meta_receptor.unique():
        # print('Looking for outliers in {}'.format(mut))
        for rate in rates_list:
            # print('Checking rate {}'.format(rate))
            rates = list(
                data.loc[data['meta_receptor'] == mut][rate])  # list() just for index reset
            cells = list(data.loc[data['meta_receptor'] == mut]['meta_file'])
            outliers_val = outliers_iqr(rates, ret='outliers')
            outliers_ind = outliers_iqr(rates, ret='outliers_indices')

            if outliers_val.size > 0:
                outliers_cel = [cells[ind] for ind in outliers_ind]
                print('{} value outliers in {} type; cells: {}, values: {} \n within: {}'.
                      format(rate, mut, outliers_cel, outliers_val, rates))
                data.loc[data['meta_file'].isin(outliers_cel), rate] = np.NAN



def statistics(data):

    all_mut_statistics = []

    for mut in data.meta_receptor.unique():
        statistics = data[data['meta_receptor'] == mut].describe()
        sem = data[data['meta_receptor'] == mut].sem(numeric_only=True, ddof=1)
        all_statistics = pd.concat([statistics, sem.to_frame().T])
        all_statistics = all_statistics.rename(index={all_statistics.index[8]: 'sem'})
        all_statistics.index = pd.MultiIndex.from_tuples([(mut, i) for i in all_statistics.index],
                                                         names=['group', 'parameter'])
        all_mut_statistics.append(all_statistics)

    return pd.concat(all_mut_statistics)

def test_statistics(data, feature_list):
    # TODO: this can be done clearer like outliers_mz

    for feature in feature_list:

        print(f'\n=== Differences analysis for "{feature}" ===\n')

        groups = [data[data['meta_receptor'] == mut][feature] for mut in data['meta_receptor'].unique()]
        cleaned_groups = [[x for x in inner_list if not pd.isna(x)] for inner_list in groups]



        print(data[['meta_receptor', feature]])
        print('')
        print(data.groupby('meta_receptor')[feature].mean())
        print('')

        # shapiro-wilk
        for mut, group in zip(data['meta_receptor'].unique(), cleaned_groups):
            if len(group) >= 3:
                s_statistic, sp_value = shapiro(group)
                if sp_value < 0.05:
                    print(""
                          "{} data not normal {}".format(mut, sp_value))
                    print(group)
                else:
                    print("{} data normal {}".format(mut, sp_value))
            else:
                print("n < 3, no Shapiro test")
        # levene
        l_statistic, lp_value = levene(*cleaned_groups)
        if lp_value < 0.05:
            print("Variances not equal {}".format(lp_value))
        else:
            print("Variances equal {}".format(lp_value))

        print('')

        if sp_value > 0.05 and lp_value > 0.05:

            f_statistic, p_value = f_oneway(*cleaned_groups)

            if p_value < 0.05:
                print("There is a significant difference between the groups (Anova f={}, p={}).".format(f_statistic, p_value))
                dropped_data = data.dropna(subset=feature)

                tuckey_results = posthoc_tukey(dropped_data, val_col=feature, group_col="meta_receptor")
                print("\nTuckey's")
                print(tuckey_results)
                dunnet_results = posthoc_dunnett(dropped_data, val_col=feature, group_col="meta_receptor", control="WT")
                print("\nDunnet's")
                print(dunnet_results)
            else:
                print("There is not a significant difference between the groups (Anova f={}, p={}).".format(f_statistic, p_value))

        else:

            h_statistic, hp_value = kruskal(*cleaned_groups)
            if hp_value < 0.05:
                print("There is a significant difference between the groups (Kruskal h={}, p={}).".format(h_statistic, hp_value))
                dropped_data = data.dropna(subset=feature)
                dunn_results = posthoc_dunn(dropped_data, val_col=feature, group_col='meta_receptor', p_adjust='holm-sidak')
                print("\nDunn's")
                print(dunn_results)
            else:
                print("There is not a significant difference between the groups (Kruskal h={}, p={}).".format(h_statistic, hp_value))




def desens_A_plot(data):

    data_desens = data.melt(id_vars='meta_receptor', value_vars=['desensitization_A%fast', 'desensitization_A%slow', 'desensitization_C%'],
                            var_name='desensitization_parameter', value_name='value')
    print(data_desens)
    sns.set_style()
    sns.set_context("paper")

    g = sns.catplot(
        data=data_desens, kind="point",
        x="desensitization_parameter", y='value', hue='meta_receptor',  estimator=np.mean, errorbar=('ci', 95), join=False,
        height=3, aspect=1.2, order=['desensitization_A%fast', 'desensitization_A%slow'],
        palette=sns.xkcd_palette(["pale red", 'windows blue', 'green']))

    g.map_dataframe(sns.swarmplot, "desensitization_parameter", 'value', alpha=.2, hue='meta_receptor',
          order=['desensitization_A%fast', 'desensitization_A%slow'],
            palette=sns.xkcd_palette(["pale red", 'windows blue', 'green']))

    g.axes[0, 0].axes.set_yticks(ticks=[0.00, 0.25, 0.5, 0.75, 1.00])
    g.ax.set_xticklabels(['A_fast', 'A_slow'])
    g.fig.set_size_inches(3.6, 3)
    g.despine(trim=True)
    g.set_axis_labels("", "")
    g._legend.set_title("mutation")

    plt.savefig('desens_A_plot' + '.png', dpi=300)


def desens_FR_plot(data):

    data_desensFR = data.melt(id_vars='meta_receptor', value_vars=['amplitudes_FR10', 'amplitudes_FR300', 'amplitudes_FR500'],
                            var_name='desensitizationFR_parameter', value_name='value')
    print(data_desensFR)
    sns.set_style()
    sns.set_context("paper")

    g = sns.catplot(
        data=data_desensFR, kind="point",
        x="desensitizationFR_parameter", y='value', hue='meta_receptor',  estimator=np.mean, errorbar=('ci', 95),  join=False,
        height=3, aspect=1.2, order=['amplitudes_FR10', 'amplitudes_FR300', 'amplitudes_FR500'],
        palette=sns.xkcd_palette(["pale red", 'windows blue', 'green']))


    g.map_dataframe(sns.swarmplot, "desensitizationFR_parameter", 'value', alpha=.2, hue='meta_receptor',
          order=['amplitudes_FR10', 'amplitudes_FR300', 'amplitudes_FR500'],
            palette=sns.xkcd_palette(["pale red", 'windows blue', 'green']))

    g.axes[0, 0].axes.set_yticks(ticks=[0.00, 0.1, 0.2, 0.3, 0.4])
    g.ax.set_xticklabels(['FR10', 'FR300', 'FR500'])
    g.fig.set_size_inches(3.6, 3)
    g.despine(trim=True)
    g.set_axis_labels("", "")
    g._legend.set_title("mutation")

    plt.savefig('desens_FR_plot' + '.png', dpi=300)


def multi_plot(data, params, params_name, labels, height, aspect, yticks):

    data_desens = data.melt(id_vars='meta_receptor', value_vars= params,
                            var_name=params_name, value_name='value')
    print(data_desens)
    sns.set_style()
    sns.set_context("paper")

    g = sns.catplot(
        data=data_desens, kind="point", dodge=0.6,
        x=params_name, y='value', hue='meta_receptor',  estimator=np.mean, errorbar=('ci', 95), join=False,
        height=height, aspect=aspect, order=params,
        palette=sns.xkcd_palette(["pale red", 'windows blue', 'green']))

    g.map_dataframe(sns.swarmplot, params_name, 'value', alpha=.2, hue='meta_receptor', dodge=True,
          order=params,
            palette=sns.xkcd_palette(["pale red", 'windows blue', 'green']))

    g.axes[0, 0].axes.set_yticks(ticks=yticks)
    g.ax.set_xticklabels(labels)
    g.fig.set_size_inches(height*aspect, height)
    g.despine(trim=True)
    g.set_axis_labels("", "")
    g._legend.set_title("mutation")

    plt.savefig(params_name + '.png', dpi=300)

def plot(feature, order, yticks=False):
    sns.set_style()
    sns.set_context("paper")

    g = sns.catplot(
        data=data, kind="point",
        x="meta_receptor", y=feature,  join=False, estimator=np.mean, errorbar=('ci', 98),
        height=2, aspect=1, order=order,
        palette=sns.xkcd_palette(["pale red", 'windows blue', 'green']
    ))


    g.map(sns.swarmplot, "meta_receptor", feature, alpha=.2, order=order,  palette=sns.xkcd_palette(["pale red", 'windows blue', 'green']))

    if yticks:
        g.axes[0, 0].axes.set_yticks(ticks=yticks)

    g.despine(trim=True)
    g.set_axis_labels("", "")
    #g.fig.suptitle(feature)
    #g.map(plt.axhline, y=1, ls='--', c='black')

    plt.savefig(feature + '.png', dpi=300)
    #plt.show()

parser = argparse.ArgumentParser()
parser.add_argument('--file_name')
args = parser.parse_args()




feature_list = ['shuts_t1', 'shuts_t2', 'shuts_t3', 'shuts_t4', 'shuts_p1', 'shuts_p2', 'shuts_p3', 'shuts_p4',
                'openings_t1',
                'rates_d', 'rates_g', 'rates_b2', 'rates_a2', 'rates_d2', 'rates_r2',
                'rates_d2p', 'rates_r2p']

'''
feature_list = ['shuts_t1', 'shuts_t2', 'shuts_t3', 'shuts_t4', 'shuts_p1', 'shuts_p2', 'shuts_p3', 'shuts_p4',
                'openings_t1', 'openings_t2',  'openings_p1', 'openings_p2',
                'rates_d', 'rates_g', 'rates_b2', 'rates_a2', 'rates_b2p', 'rates_a2p', 'rates_d2', 'rates_r2',
                'rates_d2p', 'rates_r2p']
'''

feature_list_macroscopic = ['amplitudes_RT_1090', 'amplitudes_FR10', 'amplitudes_FR300', 'amplitudes_FR500',
                            'desensitization_tau_fast', 'desensitization_tau_slow',
                            'desensitization_A%fast', 'desensitization_A%slow', 'desensitization_C%',
                            'deactivation_tau_1_comp']


data = pd.read_csv(args.file_name, header=[0, 1], sep=',')
data.columns = ['_'.join(col) for col in data.columns.values]
data.to_csv('all_data.csv')
outliers_mz(data, feature_list)
##outliers_iqrsckit(data, feature_list)
data.to_csv('no_outliers_data.csv')
statistics = statistics(data)
statistics.to_csv('statistics.csv')
test_statistics(data, feature_list)

plot('shuts_t1', ['WT', 'E153K', 'E153A'], [0.01, 0.03, 0.05, 0.07])
plot('shuts_t2', ['WT', 'E153K', 'E153A'], [0.1, 0.2, 0.3, 0.4])
plot('shuts_t3', ['WT', 'E153K', 'E153A'], [0, 1, 2, 3])
plot('shuts_t4', ['WT', 'E153K', 'E153A'], [0, 15, 30, 45, 60])

plot('openings_t1', ['WT', 'E153K', 'E153A'], [0.5, 1, 1.5, 2])
#plot('openings_t2', ['WT', 'E153K', 'E153A'], [1, 1.5, 2, 2.5, 3])


multi_plot(data, ['shuts_p1', 'shuts_p2', 'shuts_p3', 'shuts_p4'], 'shuts_distribution_parameter',
                ['P1', 'P2', 'P3', 'P4'], 3, 1.2, [0.00, 0.25, 0.5, 0.75])

#multi_plot(data, ['openings_p1', 'openings_p2'], 'openings_distribution_parameter',
#                 ['P1', 'P2'], 3, 0.8, [0.00, 0.25, 0.5, 0.75, 1])

#multi_plot(data, ['rates_d', 'rates_g', 'rates_b2', 'rates_a2', 'rates_b2p', 'rates_a2p', 'rates_d2', 'rates_r2',
#                'rates_d2p', 'rates_r2p'], 'model_rates',
#                 ['d', 'g', 'b2', 'a2', 'b2p', 'a2p', 'd', 'r', 'dp', 'rp'], 3, 2, [0,5,10,15,20])
#TODO: no explicit data passed to plot
#plot('rates_d', ['WT', 'E153K', 'E153A'], [4,5, 6, 7,])
#plot('rates_g', ['WT', 'E153K', 'E153A'], [2,5,8,11,14])
#plot('rates_b2', ['WT', 'E153K', 'E153A'], [5,15,25,35])
#plot('rates_a2', ['WT', 'E153K', 'E153A'], [0.5, 1,1.5,2, 2.5])
#plot('rates_b2p', ['WT', 'E153K', 'E153A'], [0, 3, 6, 9,])
#plot('rates_a2p', ['WT', 'E153K', 'E153A'], [0, 0.5, 1, 1.5])
#plot('rates_d2', ['WT', 'E153K', 'E153A'])
#plot('rates_r2', ['WT', 'E153K', 'E153A'])
#plot('rates_d2p', ['WT', 'E153K', 'E153A'])
#plot('rates_r2p', ['WT', 'E153K', 'E153A'])


#desens_A_plot(data)
#desens_FR_plot(data)
##plot('desensitization_A%fast', ['WT', 'E153K', 'E153A'], [.5, .75, 1])
##plot('desensitization_A%slow', ['WT', 'E153K', 'E153A'], [.0, .12, .25])
##plot('desensitization_C%', ['WT', 'E153K', 'E153A'], [.0, .12, .25])

#plot('amplitudes_RT_1090', ['WT', 'E153K', 'E153A'], [.25, .5, .75])
#plot('desensitization_tau_fast', ['WT', 'E153K', 'E153A'], [.0, 2.5, 5])
#plot('desensitization_tau_slow', ['WT', 'E153K', 'E153A'], [.0, 200, 400])
#plot('deactivation_tau_1_comp', ['WT', 'E153K', 'E153A'], [0, 250, 500])

