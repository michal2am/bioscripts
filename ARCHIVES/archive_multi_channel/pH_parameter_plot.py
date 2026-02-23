# quick script for pH paper

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px


data = pd.read_csv('amplitudes10.csv')
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

g = sns.catplot(
    data=data, kind="point",
    x="type", y="amp", hue="pH", join=False, estimator=np.mean, ci=68, dodge=True, scale=1.1,
    height=3.54, aspect=1.15,
    palette=sns.xkcd_palette(["pale red", 'windows blue']),
    order=['WT', 'E155C', 'E155S', 'E155Q', 'E155L'],
)


g.map(sns.swarmplot, "type", "amp", "pH", order=['WT', 'E155C', 'E155S', 'E155Q', 'E155L'],
      palette=sns.xkcd_palette(["pale red", 'windows blue']), alpha=0.2, marker='h')

g.axes[0, 0].axes.set_yticks(ticks=[0.5, 1, 2, 3])

g.despine(trim=True)
g.set_axis_labels("", "")
g.legend.set_title("pH")

g.map(plt.axhline, y=1, ls='--', c='black')

plt.savefig('amplitudes_ci68.png', dpi=600)
plt.show()