import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import matplotlib as mpl


class Ploter:

    def __init__(self):
        pass



    @staticmethod
    def plot_single(data, type, title):

        sns.set_context('paper')
        sns.set_style('ticks')
        font = {'family': 'serif', 'size': 12}
        mpl.rc('font', **font)

        fig, ax = plt.subplots(figsize=(3.5, 3))

        ax.plot(data)
        ax.set_xlabel(data.index.name)
        ax.set_ylabel(data.columns.values[0])

        plt.tight_layout()
        sns.despine()
        fig.savefig('test.svg', format='svg')
        plt.show()

    @staticmethod
    def plot_vertical_result(data, types, title):
        """
        plots multiple dataframes in vertical set
        :param data: dataframe top plot
        :param title: title of graph
        :return: none
        """

        dnum = len(data)
        sns.set(style='paper', font_scale=2)
        sns.set_style("ticks", {'legend.frameon': True})
        #sns.set_context("poster")

        fig = plt.figure(figsize=(2, 2), dpi=300)
        fig.canvas.set_window_title(title)
        gs = grd.GridSpec(dnum, 1, hspace=1)

        for i, (d, t) in enumerate(zip(data, types)):
            #sns.set_palette('Paired', len(d.columns))
            ax = plt.subplot(gs[i, 0])
            ax.annotate('Test',
             (d.index[2.17], d['distance [A]'][d.index[2.17]]),
             xytext=(50, 50),
             textcoords='offset points',
             arrowprops=dict(arrowstyle='simple'))
            d.plot(ax=ax, kind=t, figsize=(0.1, 0.1))
            fig.add_subplot(ax)

        #sns.despine()
        fig.savefig('test2png.png', dpi=300)
        plt.show()


