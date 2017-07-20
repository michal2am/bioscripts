import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import itertools

class Ploter:

    def __init__(self):

        # seaborn globals

        sns.set_context('paper')
        sns.set_style('ticks', {'legend.frameon': False})
        sns.set_palette("Paired", 11)

    @classmethod
    def set_rc(cls, fontsize):
        """
        sets font sizes via rcParam, use before figure generation
        :param fontsize: font size
        """

        mpl.rcParams['font.family'] = 'sans-serif'                              # global
        mpl.rcParams['font.size'] = fontsize                                    # only for relative settings

        mpl.rcParams['axes.labelsize']  = 'small'                               # element specific font size
        mpl.rcParams['axes.titlesize']  = fontsize
        mpl.rcParams['xtick.labelsize'] = fontsize
        mpl.rcParams['ytick.labelsize'] = fontsize
        mpl.rcParams['legend.fontsize'] = 'x-small'

    # figure handling

    @classmethod
    def prepare_single_figure(cls, title, size):
        fig, ax = plt.subplots(figsize=size)
        fig.canvas.set_window_title(title)
        return fig, ax

    @classmethod
    def prepare_multiple_cx_figure(cls, title, size, length):
        fig, axes = plt.subplots(figsize=size, nrows=length, sharex=True)
        fig.canvas.set_window_title(title)
        return fig, axes

    @classmethod
    def save_figure(cls, figure, title, rect):
        plt.tight_layout(rect=rect)
        sns.despine()
        figure.savefig(title+'.svg', format='svg')
        plt.show()

    # axes handling

    @classmethod
    def set_axes(cls, ax, **kwargs):
        """
        :param ax:
        :param kwargs:
        :return:
        """

        if 'xtformat' in kwargs: ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(kwargs['xtformat']))
        if 'ytformat' in kwargs: ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(kwargs['ytformat']))

        if 'xtvals' in kwargs: ax.xaxis.set_ticks(kwargs['xtvals'])
        if 'xtlabs' in kwargs: ax.xaxis.set_ticklabels(kwargs['xtlabs'])
        if 'ytvals' in kwargs: ax.yaxis.set_ticks(kwargs['ytvals'])
        if 'ytlabs' in kwargs: ax.yaxis.set_ticklabels(kwargs['ytlabs'])

        if 'xlim' in kwargs: ax.set_xlim(kwargs['xlim'])
        if 'ylim' in kwargs: ax.set_ylim(kwargs['ylim'])

        #ax.locator_params(axis='x', nbins=5)
        #ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))

    @classmethod
    def set_legend(cls, ax, **kwargs):
        """
        sets legend parameters
        :param ax: axes
        """

        leg = ax.legend(**kwargs)
        for legobj in leg.legendHandles:
            legobj.set_linewidth(2.0)

    @classmethod
    def set_annotation_ranges(cls, ax, data, **kwargs):

        if set(kwargs['features']).intersection(set(list(data.columns))):

            i = 0.08
            annots = []
            for feature in kwargs['features']:

                data[feature] = data[feature].replace('-', np.nan)
                annots.append([i for i in data[feature].unique() if not pd.isnull(i)])
                data = pd.concat([data, pd.get_dummies(data[feature]).replace(0, np.nan) * -i ])
                i += 0.002

            return data, annots

        else:
            return data, None

    @classmethod
    def set_annotation_legend(cls, ax, annots, kwargs):
        """
        
        :param ax: 
        :param annots: 
        :param kwargs: 
        :return: 
        """
        annots = set(list(itertools.chain(*annots)))
        handles, labels = ax.get_legend_handles_labels()
        labels = set(labels) - annots
        kwargs['legend_style']['labels'] = labels

    @staticmethod
    def plot_single(style, data, title, size, kind='line', font=18, **kwargs):
        """
        plot single plot of one dataseries
        :param style: 'single' or 'multi'
        :param data: pandas dataframe
        :param title: string
        :param size: tuple, figure size in inches
        :param kind:
        :param font: int, font size
        """

        # figure and axes initialization
        Ploter.set_rc(font)
        fig, ax = Ploter.prepare_single_figure(title, size)

        # for single dataframe
        if style == 'single':

            data.plot(ax=ax, kind=kind, **kwargs.get('lines_style', {}))
            ax.set_xlabel(data.index.name)
            if len(data.columns.values) == 1:
                ax.set_ylabel(data.columns.values[0])
            else:
                ax.set_ylabel(kwargs['y_label'])
                Ploter.set_legend(ax, **kwargs.get('legend_style', {}))

        # for list of dataframes
        if style == 'multi':
            for series in data:

                if 'annotation' in kwargs:
                    series, annots = Ploter.set_annotation_ranges(ax, series, **kwargs.get('annotation'))

                series.plot(ax=ax, kind=kind, **kwargs.get('lines_style', {}))

            if 'annotation' in kwargs: Ploter.set_annotation_legend(ax, annots, kwargs)

            ax.set_xlabel(kwargs['x_label'])
            ax.set_ylabel(kwargs['y_label'])

            Ploter.set_legend(ax, **kwargs.get('legend_style', {}))

        Ploter.set_axes(ax, **kwargs.get('axes_style', {}))
        Ploter.save_figure(fig, title, kwargs.get('rect', (0, 0, 1, 1)))
