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
    data=data, kind="box",
    x="type", y="amp", hue="pH",
    ci="sd", height=3.54, aspect=1.15,
    order=['wt', 'cys', 'ser', 'gln', 'leu']
)

#g.map(sns.swarmplot, "type", "amp", "pH", order=['wt', 'cys', 'ser', 'gln', 'leu'])


g.despine()
g.set_axis_labels("", "")
g.legend.set_title("")

g.map(plt.axhline, y=1, ls='--', c='black')

plt.savefig('amplitudes.png')
plt.show()