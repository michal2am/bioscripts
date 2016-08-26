import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd


class Ploter:

    def __init__(self):
        pass

    @staticmethod
    def plot_vertical_result(data, types, title):
        """
        plots multiple dataframes in vertical set
        :param data: dataframe top plot
        :param title: title of graph
        :return: none
        """

        dnum = len(data)

        sns.set_style("ticks", {'legend.frameon': True})
        sns.set_context("poster")

        fig = plt.figure()
        fig.canvas.set_window_title(title)
        gs = grd.GridSpec(dnum, 1, hspace=0.5)

        for i, (d, t) in enumerate(zip(data, types)):
            sns.set_palette('Paired', len(d.columns))
            ax = plt.subplot(gs[i, 0])
            d.plot(ax=ax, kind=t)
            fig.add_subplot(ax)

        sns.despine()
        plt.show()


