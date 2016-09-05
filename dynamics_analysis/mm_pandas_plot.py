import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import matplotlib as mpl


class Ploter:

    def __init__(self):
        pass

    @classmethod
    def extra_annotation(cls, ax, text, xypoint, xytext):
        # ax.annotate(text, xypoint, xycoords='data', xytext=xytext, textcoords='data', arrowprops=dict(arrowstyle='simple', color='black'))
        ax.annotate(text, xypoint, arrowprops=dict(arrowstyle='simple', color='black'))

    @staticmethod
    def plot_single(data, title, **kwargs):
        """
        plot single plot of one dataseries
        :param data: pandas dataframe
        :param title: string
        :return:
        """

        # seaborn globals
        sns.set_context('paper')
        sns.set_style('ticks')
        font = {'family': 'serif', 'size': 12}
        mpl.rc('font', **font)

        # figure and axes initialization
        fig, ax = plt.subplots(figsize=(4, 3))
        fig.canvas.set_window_title(title)

        # data and labels selection
        ax.plot(data, color='black', linewidth=0.5)
        ax.set_xlabel(data.index.name)
        ax.set_ylabel(data.columns.values[0])

        # annotations
        points = (2.17, 5.76)
        annotations = (annotation for annotation in kwargs['annotations'] if 'annotations' in kwargs.keys())
        for annotation in annotations:
            Ploter.extra_annotation(ax, 'protein constraints removed', (points[0], data.loc[points[0], data.columns.values[0]]), (0.1, 0.7))
            Ploter.extra_annotation(ax, 'ion constraints removed', (points[1], data.loc[points[1], data.columns.values[0]]), (0.1, 0.85))

        # ticks formatting
        ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
        ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
        ax.xaxis.set_ticks(np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], num=5))
        ax.yaxis.set_ticks(np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], num=5))

        # font formatting
        font = {'family': 'serif'}
        mpl.rc('font', **font)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(12)

        # post plot mods and saving
        plt.tight_layout()
        sns.despine()
        fig.savefig(title+'.svg', format='svg')
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


