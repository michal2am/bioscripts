import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px

data = pd.read_csv('multiplot.csv')
data_long = pd.read_csv('multiplot_long.csv')
data_long_selected = pd.read_csv('multiplot_long_selected.csv')

print(data)

sns.set_style()
sns.set_context("paper")

#g = sns.catplot(kind='point', data=data, x='pH', y='amp',
#                hue='type', dodge=True, legend_out=True)

# g.set_titles("{col_name} {row_name}")
#g.set_titles("{col_name}")

#g.set_axis_labels("", "")

#plt.tight_layout()
#plt.savefig('seaborn_plot_{}_lignad.png'.format(chain))
#plt.show()

g = sns.relplot(
    data=data,
    x="Tyr97", y="Ser156", hue="type", #col='pose'
    #col_order=['acidic', 'basic'],
    #row='type', col='pH_range',
    #height=2, aspect=1,
    #palette=sns.xkcd_palette(["pale red", "greyish", 'windows blue'])
)

print(data)
data_agg = data.groupby('type').mean()

g = sns.relplot(
    data=data_agg,
    x="Tyr97", y="Tyr157", kind='line', hue="type"#, col='pose'
    #col_order=['acidic', 'basic'],
    #row='type', col='pH_range',
    #height=2, aspect=1,
    #palette=sns.xkcd_palette(["pale red", "greyish", 'windows blue'])
)
#g = sns.PairGrid(data, hue="type")
#g.map_diag(sns.histplot)
#g.map_offdiag(sns.scatterplot)

#g.map_offdiag(sns.scatterplot, size=data[data['pose']])

'''
g = sns.catplot(
    data=data_long, kind='box',
    x="type", y="distance", hue='residue',
    #col_order=['acidic', 'basic'],
    #row='type', col='pH_range',
    #height=2, aspect=1,
    #palette=sns.xkcd_palette(["pale red", "greyish", 'windows blue'])
    order=['WT', 'E155C', 'E155S', 'E155Q', 'E155L', 'E155A', 'E155N'],

)
'''

#PAPER

g = sns.catplot(
    data=data_long_selected, kind='box',
    x="type", y="distance", hue='residue',
    ci="sd", height=3.54, aspect=1.15,
    #col_order=['acidic', 'basic'],
    #row='type', col='pH_range',
    #height=2, aspect=1,
    #palette=sns.xkcd_palette(["pale red", "greyish", 'windows blue'])
    palette=sns.xkcd_palette(['windows blue','dusty orange', 'dark pastel green']),

    order=['WT', 'E155C', 'E155S', 'E155Q', 'E155L'],

)

#PAPER TOO
#g.axes[0, 0].axes.set_yticks(ticks=[2.5, 3, 3.5, 4.5])
g.axes[0, 0].yaxis.set_major_locator(plt.MaxNLocator(5))

g.despine(trim=True)
g.set_axis_labels("", "")
g.legend.set_title("distance to")

plt.savefig('distances.png', dpi=300)
plt.show()


g = sns.displot(
    data=data_long_selected, #kind='box',
    x="distance", row='residue', hue='type', kind='kde',
    #col_order=['acidic', 'basic'],
    #row='type', col='pH_range',
    #height=2, aspect=1,
    #palette=sns.xkcd_palette(["pale red", "greyish", 'windows blue'])
    #order=['WT', 'E155C', 'E155S', 'E155Q', 'E155L', 'E155A', 'E155N'],

)

g = sns.displot(
    data=data, #kind='box',
    x="Tyr97", y="Tyr157", hue='type', kind='kde',
    #col_order=['acidic', 'basic'],
    #row='type', col='pH_range',
    #height=2, aspect=1,
    #palette=sns.xkcd_palette(["pale red", "greyish", 'windows blue'])
    #order=['WT', 'E155C', 'E155S', 'E155Q', 'E155L', 'E155A', 'E155N'],

)

#g.map(sns.swarmplot, "type", "amp", "pH", order=['wt', 'cys', 'ser', 'gln', 'leu'])


#g.despine()
#g.set_axis_labels("", "")
#g.legend.set_title("")

#g.map(plt.axhline, y=1, ls='--', c='black')

#plt.savefig('distances.png')
plt.show()