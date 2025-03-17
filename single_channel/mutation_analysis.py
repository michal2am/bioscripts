import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import f_oneway
from scipy.stats import shapiro
from scipy.stats import levene
from scikit_posthocs import outliers_iqr
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import argparse



sns.set_style()
sns.set_context("talk")
sns.set_palette('muted')


def outliers_iqrsckit(data, rates_list):
    for mut in data.meta_residue_mut.unique():
        # print('Looking for outliers in {}'.format(mut))
        for rate in rates_list:
            # print('Checking rate {}'.format(rate))
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


def statistics(data):

    all_mut_statistics = []

    for mut in data.meta_residue_mut.unique():
        statistics = data[data['meta_residue_mut'] == mut].describe()
        sem = data[data['meta_residue_mut'] == mut].sem(numeric_only=True, ddof=1)
        all_statistics = pd.concat([statistics, sem.to_frame().T])
        all_statistics = all_statistics.rename(index={all_statistics.index[8]: 'sem'})
        all_statistics.index = pd.MultiIndex.from_tuples([(mut, i) for i in all_statistics.index],
                                                         names=['group', 'parameter'])
        all_mut_statistics.append(all_statistics)

    return pd.concat(all_mut_statistics)

def test_statistics(data, feature_list):
    with open('test_statistics.log', 'w') as f:
        for feature in feature_list:
            groups = [data[data['meta_residue_mut'] == mut][feature] for mut in data['meta_residue_mut'].unique()]
            cleaned_groups = [[x for x in inner_list if not pd.isna(x)] for inner_list in groups]
            print(feature)  # , file=f)
            for group in cleaned_groups:
                if len(group) >= 3:
                    s_statistic, sp_value = shapiro(group)
                    if sp_value < 0.05:
                        print("Data not normal {}".format(sp_value))
                        print(group)
                    else:
                        print("Data normal {}".format(sp_value))
                    l_statistic, lp_value =
                else:
                    print("n < 3, no Shapiro test")
            f_statistic, p_value = f_oneway(*cleaned_groups)
            print(f_statistic, p_value)  # , file=f)
            if p_value < 0.05:
                print("There is a significant difference between the groups.")  # , file=f)
                dropped_data = data.dropna(subset=feature)
                tukey = pairwise_tukeyhsd(endog=dropped_data[feature], groups=dropped_data['meta_residue_mut'],
                                          alpha=0.05)
                print(tukey.summary())  # , file=f)
            else:
                print("There is not a significant difference between the groups.")  # , file=f)

parser = argparse.ArgumentParser()
parser.add_argument('--file_name')
args = parser.parse_args()

feature_list = ['shuts_t1', 'shuts_t2', 'shuts_t3', 'shuts_t4', 'shuts_p1', 'shuts_p2', 'shuts_p3', 'shuts_p4',
                'openings_t1', 'openings_t2',  'openings_p1', 'openings_p2',
                'rates_d', 'rates_g', 'rates_b2', 'rates_a2', 'rates_b2p', 'rates_a2p', 'rates_d2', 'rates_r2',
                'rates_d2p', 'rates_r2p']


data = pd.read_csv(args.file_name, header=[0, 1], sep=';')
data.columns = ['_'.join(col) for col in data.columns.values]
data.to_csv('all_data.csv')
outliers_iqrsckit(data, feature_list)
data.to_csv('no_outliers_data.csv')
statistics = statistics(data)
statistics.to_csv('statistics.csv')
test_statistics(data, feature_list)

