import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px


data = pd.read_csv('amplitudes1030.csv')
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
'''
g = sns.catplot(
    data=data, kind="point",
    x="type", y="amp", hue="pH", join=False, estimator=np.mean, ci=95, dodge=True,
    height=3.54, aspect=1.15,
    palette=sns.xkcd_palette(["pale red", 'windows blue']),
    order=['WT30', 'WT10', 'E155C', 'E155S', 'E155Q', 'E155L'],
)

g.map(sns.swarmplot, "type", "amp", "pH", order=['WT30', 'WT10', 'E155C', 'E155S', 'E155Q', 'E155L'],
      palette=sns.xkcd_palette(["pale red", 'windows blue']), alpha=0.5, marker='h')
'''
g = sns.catplot(
    data=data, kind="point",
    x="pH", y="amp", hue="type", join=True, estimator=np.mean, ci=95, dodge=False,
    height=3.54, aspect=1.15,
    #palette=sns.xkcd_palette(["pale red", 'windows blue']),
    hue_order=['WT10', 'E155C', 'E155S', 'E155Q', 'E155L']
)

#g.map(sns.swarmplot, "pH", "amp", "type", hue_order=['WT10', 'E155C', 'E155S', 'E155Q', 'E155L'],
#       alpha=0.5, marker='h')




g.axes[0, 0].axes.set_yticks(ticks=[0.5, 1, 2, 3])

g.despine(trim=True)
g.set_axis_labels("", "")
g.legend.set_title("")

g.map(plt.axhline, y=1, ls='--', c='black')

plt.savefig('amplitudes.png', dpi=600)
plt.show()