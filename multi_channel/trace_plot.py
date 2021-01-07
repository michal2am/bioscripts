import pyabf
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import pandas as pd
import seaborn as sns

data = pd.read_csv('trace_config.csv')
print(data)
read_data = pd.DataFrame()

for idx, abf_row in data.iterrows():

    print(abf_row['abf'])

    read_abf = pyabf.ABF(abf_row['abf'])

    new_cell = pd.DataFrame()

    for sweep_no, sweep in enumerate(range(read_abf.sweepCount)):
        read_abf.setSweep(sweep)
        # print(sweep_no)
        new_cell['sweep_' + str(sweep_no)] = read_abf.sweepY[0:1010]
        t = read_abf.sweepX[0:1010]
        # print(read_abf.sweepY)
        #plt.plot(read_abf.sweepX, read_abf.sweepY)
        #plt.show()

    new_cell['sweep_mean'] = new_cell.mean(axis=1)
    new_cell['t'] = t *1000
    new_cell['t'] += abf_row['shift']
    new_cell['type'] = abf_row['type']
    new_cell['pH'] = abf_row['pH']
    new_cell['pH_range'] = abf_row['pH_range']
    # print(new_cell)
    read_data = pd.concat([read_data, new_cell])

print(read_data)



for receptor in ['WT', 'E155C', 'E155S', 'E155Q', 'E155L']:

    receptor_data = read_data[read_data.loc[:, 'type'] == receptor]

    sns.set_style()
    sns.set_context("paper")

    g = sns.relplot(
        data=receptor_data, kind="line",
        x="t", y="sweep_mean",  hue="pH",
        col_order=['acidic', 'basic'],
        row='type', col='pH_range',
        height=2, aspect=1,
        palette=sns.xkcd_palette(["pale red", "greyish", 'windows blue'])
    )


    #g.map(sns.swarmplot, "type", "amp", "pH", order=['wt', 'cys', 'ser', 'gln', 'leu'])
    g.set_titles(col_template="{col_name}", row_template="{row_name}")

    #g.despine(top=False, bottom=True, trim=True)
    g.despine()

    g.set_axis_labels("", "")
    g.legend.set_title("pH")
    #g.legend.remove()

    # tick number setters
    g.axes[0, 0].xaxis.set_major_locator(plt.MaxNLocator(5))
    #g.axes[0, 0].yaxis.set_major_locator(plt.MaxNLocator(4))

    # fixed axes limits and ticks
    #g.axes[0, 0].axes.set_xlim(0, None)
    g.axes[0, 0].axes.set_ylim(None, 25  )
    g.axes[0, 0].axes.set_yticks(ticks=[-800, -600, -400, -200, 0])

    plt.tight_layout()
    plt.savefig('traces_{}.png'.format(receptor), dpi=600)
    plt.show()