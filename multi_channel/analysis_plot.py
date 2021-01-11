import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px

data = pd.read_csv('smarts_mm.csv')
print(data)

data_long = pd.melt(data, id_vars=['lp', 'type', 'configuration', 'file', 'pulse_length', 'sweep', 'concentration'])
print(data_long)


sns.set_style()
sns.set_context('talk')



def des_joint_plot():

    for des in [['des_Af', 'des_tf', 'des_f'], ['des_As', 'des_ts', 'des_s']]:

        g = sns.displot(data=data, kind='kde',
                        x=des[0], y=des[1], hue='type',
                        levels=1,
                        hue_order=['WT', 'L300V', 'L296V', 'L300V+L296V'],
                        height=4, aspect=1.75,
                        palette=sns.color_palette('deep', n_colors=4),
                        )
        g.map(sns.scatterplot, des[0], des[1], 'type',
              hue_order=['WT', 'L300V', 'L296V', 'L300V+L296V'],
              palette=sns.color_palette('deep', n_colors=4),
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
        order=['WT', 'L300V', 'L296V', 'L300V+L296V'],
        hue_order=['WT', 'L300V', 'L296V', 'L300V+L296V'],
        height=4, aspect=1.75,
        palette=sns.color_palette('deep'),
    )

    g.map(sns.swarmplot, "type", param, 'type',
          order=['WT', 'L300V', 'L296V', 'L300V+L296V'],
          hue_order=['WT', 'L300V', 'L296V', 'L300V+L296V'],
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
        hue_order=['WT', 'L300V', 'L296V', 'L300V+L296V'],
        height=4, aspect=1.75,
        palette=sns.color_palette('deep'),
    )

    g.map(sns.swarmplot, 'variable', 'value', 'type',
          order=['FR10', 'FR300', 'FR500'],
          hue_order=['WT', 'L300V', 'L296V', 'L300V+L296V'],
          palette=sns.color_palette('deep'),
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
        hue_order=['WT', 'L300V', 'L296V', 'L300V+L296V'],
        height=4, aspect=1.75,
        palette=sns.color_palette('deep'),
    )

    g.map(sns.swarmplot, 'variable', 'value', 'type',
          order=['des_Af', 'des_As', 'des_AC'],
          hue_order=['WT', 'L300V', 'L296V', 'L300V+L296V'],
          palette=sns.color_palette('deep'),
          alpha=0.5, marker='h')

    g.axes[0, 0].axes.set_yticks(ticks=[0, .25, .5, .75, 1])
    #g.axes[0, 0].yaxis.set_major_locator(plt.MaxNLocator(5))

    g.despine(trim=True)
    g.set_axis_labels("", "")
    g.legend.set_title("receptor type")

    plt.savefig('des_a_all.png', dpi=600)
    plt.show()






'''
# jointy dla des

'''


'''
# joint przez pairgrida
g = sns.PairGrid(data=data, vars=["des_Af","des_tf", 'des_As', 'des_ts'], hue='type',)
g.map_upper(sns.scatterplot)
g.map_lower(sns.kdeplot)
g.map_diag(sns.kdeplot, lw=1, legend=False)
'''

'''
def des_joint_plot_full():

    g = sns.jointplot(
        data=data, x="des_Af", y="des_tf", hue='type',
                      hue_order=['WT', 'L300V', 'L296V', 'L300V+L296V'],
                      kind='scatter'
                      )
    g.plot_joint(sns.kdeplot, levels=1)
    # g.plot_marginals(sns.histplot)
    # sns.despine(trim=True)
    plt.show()

    g = sns.jointplot(data=data, x="des_As", y="des_ts", hue='type',
                      hue_order=['WT', 'L300V', 'L296V', 'L300V+L296V'],
                      kind='scatter'
                      )
    g.plot_joint(sns.kdeplot, levels=1)
    # g.plot_marginals(sns.histplot)
    # sns.despine(trim=True)
    plt.show()
'''

# pairplot na wszystko
# sns.pairplot(data, hue='mutant')
# plt.tight_layout()
# plt.show()


#for param in ['rt_10/90', 'FR10', 'FR300', 'FR500', 'des_ts', 'des_tf', 'des_As', 'des_Af', 'des_AC', 'dea_ts', 'dea_tf', 'dea_As', 'dea_Af', 'dea_tm']:
#    interval_plot(param)

interval_plot('rt_1090')
interval_plot('dea_tm')
interval_plot('des_tf')
interval_plot('des_ts')

fr_plot()

des_joint_plot()
des_a_plot()

