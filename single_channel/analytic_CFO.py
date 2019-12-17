import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import itertools as itr
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


class CFO:

    def __init__(self, rates):

        self.delta, self.gamma, self.beta, self.alpha = rates

        self.ts1, self.ts2, self.ts3, self.p1, self.p2, self.p3 = self.calculate_equiprobs()
        self.r, self.r_app, self.l1, self.l2, self.t1, self.t2, self.a1, self.a2, self.mb, \
        self.a1_app, self.a2_app, self.mb_app, \
        self.mopb, self.mopb_app, self.mshb, self.mshb_app = self.calculate_bursts()
        self.l3, self.l4, self.t3, self.t4, self.a3, self.a4 = self.calculate_shuts()

        to_frame = {'delta': [self.delta], 'gamma': [self.gamma], 'beta': [self.beta], 'alpha': [self.alpha],
                    'ts1': [self.ts1], 'ts2': [self.ts2], 'ts3': [self.ts3],  # states' life times
                    'p1': [self.p1], 'p2': [self.p2], 'p3': [self.p3],  # states' equilibrium occupancies
                    'r': [self.r], 'r_app': [self.r_app],  # real/apparent openings per burst
                    't1': [self.t1], 't2': [self.t2], 'a1': [self.a1], 'a2': [self.a2], 'mb': [self.mb],  # real burst time distribution
                    'a1_app': [self.a1_app], 'a2_app': [self.a2_app], 'mb_app': [self.mb_app],  # apparent burst time distribution
                    'mopb': [self.mopb], 'mopb_app': [self.mopb_app], 'mshb': [self.mshb], 'mshb_app': [self.mshb_app],  # real/apparent mean open and shut times per birst
                    't3': [self.t3], 't4': [self.t4], 'a3': [self.a3], 'a4': [self.a4],  # shut time distribution
                    }

        self.frame = pd.DataFrame(to_frame)

    def calculate_equiprobs(self):

        # equilibrium occupancies:
        p3 = (((self.delta/self.gamma)*(1+self.beta/self.alpha)+1)**-1)*100
        p2 = ((self.delta/self.gamma)*p3/100)*100
        p1 = ((self.delta/self.gamma)*(self.beta/self.alpha)*p3/100)*100

        # mean state life time
        ts1 = (1/self.alpha)*1000
        ts2 = (1/(self.beta + self.gamma))*1000
        ts3 = (1/self.delta)*1000

        return ts1, ts2, ts3, p1, p2, p3

    def calculate_bursts(self):

        # mean and mean apparent openings per burst
        r = self.beta/self.gamma
        r_app =1 + self.beta/self.gamma

        # burst length distribution
        b = self.alpha + self.beta + self.gamma
        c = self.alpha * self.gamma

        l1 = 0.5*(b - np.sqrt(b**2 - 4*c))
        l2 = 0.5*(b + np.sqrt(b**2 - 4*c))

        t1 = l1**-1*1000
        t2 = l2**-1*1000

        const = self.gamma/(l2 - l1)

        a1 = const * (self.alpha - l1)
        a2 = const * (l2 - self.alpha)

        mb = (self.alpha + self.beta)/(self.alpha*self.gamma)*1000

        const_app = l1*l2/((l2-l1)*(self.beta+self.gamma))

        a1_app = const_app * (self.beta + self.gamma - l1)
        a2_app = const_app * (l2 - self.beta - self.gamma)

        mb_app = ((1+self.beta/self.gamma)*1/self.alpha+self.beta/self.gamma*(1/(self.beta+self.gamma)))*1000

        # mean open and shut time per burst

        mopb = (self.beta/(self.gamma*self.alpha))*1000
        mopb_app = ((self.beta+self.gamma)/(self.gamma*self.alpha))*1000

        mshb = (1/self.gamma)*1000
        mshb_app = (self.beta/(self.gamma*(self.beta+self.gamma)))*1000

        return r, r_app, l1, l2, t1, t2, a1, a2, mb, a1_app, a2_app, mb_app, mopb, mopb_app, mshb, mshb_app

    def calculate_shuts(self):

        # shut time distributions
        b = self.beta + self.delta + self.gamma
        c = self.beta * self.delta

        l3 = 0.5*(b - np.sqrt(b**2 - 4*c))
        l4 = 0.5*(b + np.sqrt(b**2 - 4*c))

        t3 = l3**-1*1000
        t4 = l4**-1*1000

        const = self.beta/(l4 - l3)

        a3 = const*(self.delta - l3)
        a4 = const*(l4 - self.delta)

        return l3, l4, t3, t4, a3, a4


def rate_range(rate):
    #step = int(rate * 0.1)
    #return list(range(rate - 9 * step, rate + 10 * step, step))
    #return list(np.linspace(100, 2500, 20))
    return [int(rate) for rate in list(np.logspace(2.5, 4, num=15))]

def rate_combinations(rates):
    return list(itr.product(*rates))


dc = CFO([26, 1e4, 1.9e4, 1e3])
print(dc.frame)
print('DC sanity: 58.7 ms, 34.5 us, 5.879, 18994')


def generate():

    delta, gamma, beta, alpha = 1500, 2000, 8200, 700
    deltas, gammas, betas, alphas = rate_range(delta), rate_range(gamma), rate_range(beta), rate_range(alpha)
    rate_sets = rate_combinations([deltas, gammas, betas, alphas])                                                          # 14641 combinations

    tests = []

    for rate_set in rate_sets:

        cfo_mechanism = CFO(rate_set)
        tests.append(cfo_mechanism.frame)

    all_tests = pd.concat(tests, axis=0)
    all_tests.reset_index(inplace=True, drop=True)
    print(all_tests)

    results_long = all_tests.melt(id_vars=['delta', 'gamma', 'beta', 'alpha']).copy()

    all_tests.to_csv('all_test_log_mini.csv')


delta, gamma, beta, alpha = 1500, 2000, 8200, 700

#generate()
all_tests = pd.read_csv('all_test_log_mini.csv', index_col=0)


'''
fix_not_delta = all_tests[((all_tests['gamma'] == gamma) & (all_tests['beta'] == beta) & (all_tests['alpha'] == alpha))]
fix_not_gamma = all_tests[((all_tests['delta'] == delta) & (all_tests['beta'] == beta) & (all_tests['alpha'] == alpha))]
fix_not_beta = all_tests[((all_tests['delta'] == delta) & (all_tests['gamma'] == gamma) & (all_tests['alpha'] == alpha))]
fix_not_alpha = all_tests[((all_tests['delta'] == delta) & (all_tests['beta'] == beta) & (all_tests['beta'] == beta))]


fix_beta_alpha = all_tests[((all_tests['beta'] == beta) & (all_tests['alpha'] == alpha))]
fix_delta_gamma = all_tests[((all_tests['gamma'] == gamma) & (all_tests['delta'] == delta))]

fix_gamma_alpha = all_tests[((all_tests['gamma'] == gamma) & (all_tests['alpha'] == alpha))]
fix_gamma_beta = all_tests[((all_tests['gamma'] == gamma) & (all_tests['beta'] == beta))]

#fix_delta_alpha = all_tests[((all_tests['delta'] == delta) & (all_tests['alpha'] == alpha))]
fix_delta_beta = all_tests[((all_tests['delta'] == delta) & (all_tests['alpha'] == alpha))]

sns.set_context("talk")

'''


def correlation_plot():

    correlations_pearson = all_tests.corr()
    annots_pearson = correlations_pearson.iloc[4:, 0:4].round(decimals=1)
    annots_pearson.replace(to_replace=[0.0, -0.0], value='', inplace=True)

    correlations_kendall = all_tests.corr(method='kendall')
    annots_kendall = correlations_kendall.iloc[4:, 0:4].round(decimals=1)
    annots_kendall.replace(to_replace=[0.0, -0.0], value='', inplace=True)

    correlations_spearman = all_tests.corr(method='spearman')
    annots_spearman = correlations_spearman.iloc[4:, 0:4].round(decimals=1)
    annots_spearman.replace(to_replace=[0.0, -0.0], value='', inplace=True)

    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=[15, 5])
    sns.heatmap(correlations_pearson.iloc[4:, 0:4], yticklabels=1, cmap="RdBu_r", center=0, annot=annots_pearson, fmt='', ax=axs[0])
    sns.heatmap(correlations_kendall.iloc[4:, 0:4], yticklabels=1, cmap="RdBu_r", center=0, annot=annots_kendall, fmt='', ax=axs[1])
    sns.heatmap(correlations_spearman.iloc[4:, 0:4], yticklabels=1, cmap="RdBu_r", center=0, annot=annots_spearman, fmt='', ax=axs[2])

# correlation_plot()


def two_rate_heatmap(rates, fix_rates, property):

    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=[5, 5])

    fixed = all_tests[((all_tests[fix_rates[0]] == 316) & (all_tests[fix_rates[1]] == 316))]
    for_heat = fixed[[rates[0], rates[1], property]].pivot(index=rates[0], columns=rates[1], values=property)
    for_heat = for_heat.reindex(index=for_heat.index[::-1])

    ax.set_title(property)
    sns.heatmap(for_heat, annot=False, xticklabels=2, yticklabels=2, cbar=True, square=True, ax=ax)


# two_rate_heatmap(['gamma', 'beta'], ['alpha', 'delta'], 'r')
#two_rate_heatmap(['gamma', 'beta'], ['alpha', 'delta'], 'r_app')
#two_rate_heatmap(['gamma', 'beta'], ['alpha', 'delta'], 'mshb')
#two_rate_heatmap(['gamma', 'beta'], ['alpha', 'delta'], 'mshb_app')




# TODO: add rates parameter and omit parameter, basically rewrite for 3-rate dependancy option

def four_rate_plot(property):

    fig, axs = plt.subplots(ncols=2, nrows=2, figsize=[20, 20], sharey=True)
    flat_axs = [ax for col in axs for ax in col]

    for rate, ax in zip(['delta', 'gamma', 'beta', 'alpha'], flat_axs):

        sns.boxenplot(x=rate, y=property, data=all_tests, ax=ax, color='grey')
        sns.stripplot(x=rate, y=property, data=all_tests, ax=ax, marker='.', size=3, jitter=True, color='black')

        sns.despine()


def one_rate_plots(property):

    rates = ['delta', 'gamma', 'beta', 'alpha']

    for rate in rates:

        fig, axs = plt.subplots(ncols=1, nrows=3, figsize=[20, 20])
        hue_rates = rates.copy()
        hue_rates.remove(rate)

        for hue_rate, ax in zip(hue_rates, axs):
            sns.stripplot(x=rate, y=property, data=all_tests, ax=ax, marker='.', color='red', size=3,
                                hue=hue_rate, jitter=True, dodge=True)





#four_rate_plot('a4')
#four_rate_plot('a3')
#four_rate_plot('t4')
#four_rate_plot('t3')


#one_rate_plots('r')

sns.scatterplot(x='delta', y='gamma', data=all_tests, style='beta', hue='a4')

'''
def plot_paired(models, factor):

    g = sns.PairGrid(models, vars=['t3', 't4', 'a3', 'a4'], hue=factor)#, palette=sns.color_palette("coolwarm", 20))
    g = g.map_offdiag(plt.scatter)
    g = g.add_legend()
    sns.despine(trim=True)
plot_paired(selected_gamma_delta, 'delta')
'''

#sns.set_style("white")
#sns.set_context("talk")

sns.despine()
plt.tight_layout()
plt.show()

