import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import matplotlib as mpl


class Ploter:

    def __init__(self):

        # seaborn globals

        sns.set_context('paper')
        sns.set_style('ticks', {'legend.frameon': True})
        font = {'family': 'serif', 'size': 12}
        mpl.rc('font', **font)

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
    def save_figure(cls, figure, title):
        plt.tight_layout()
        sns.despine()
        figure.savefig(title+'.svg', format='svg')
        plt.show()

    # axes handling

    @classmethod
    def set_ticks(cls, ax, freq):
        ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
        ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
        ax.xaxis.set_ticks(np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], num=freq))
        ax.yaxis.set_ticks(np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], num=freq))

        ax.xaxis.set_ticks(np.linspace(0, 20, num=freq))
        ax.set_xlim([-5, 30])

    @classmethod
    def set_fonts(cls, ax, size):
        font = {'family': 'serif'}
        mpl.rc('font', **font)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(size)

    @classmethod
    def set_legend(cls, ax, position, fontsize, cols=1):
        leg = ax.legend(fontsize=fontsize, loc=position, ncol=cols)
        for legobj in leg.legendHandles:
            legobj.set_linewidth(2.0)

    @classmethod
    def set_annotation(cls, ax, text, xypoint): #, xytext):
        # TODO: xytext!
        # ax.annotate(text, xypoint, xycoords='data', xytext=xytext, textcoords='data', arrowprops=dict(arrowstyle='simple', color='black'))
        ax.annotate(text, xypoint, xycoords='data', textcoords='data', arrowprops=dict(arrowstyle='simple', color='black'))

    @staticmethod
    def plot_single(style, data, title, size, font=12, ticks=5, **kwargs):
        """
        plot single plot of one dataseries
        :param style: 'single' or 'multi-indepx'
        :param data: pandas dataframe
        :param title: string
        :param size: tuple, figure size in inches
        :param font: int, font size
        :param ticks: int, tick plot number
        """

        # kwargs clean up
        kwargs = dict((k, v) for k, v in kwargs.items() if kwargs[k] is not None)

        # figure and axes initialization
        fig, ax = Ploter.prepare_single_figure(title, size)

        # data and labels selection

        # for dataframe with one column
        if style == 'single':
            ax.plot(data, 'o', linewidth=1.5)
            ax.set_xlabel(data.index.name)
            ax.set_ylabel(data.columns.values[0])

        # for dataframe with many columns
        if style == 'multi':
            data.plot(ax=ax, linewidth=0.25)                                    # correct label and legend handling
            ax.set_xlabel(data.index.name)
            ax.set_ylabel(kwargs['y_label'])
            Ploter.set_legend(ax, 4, font, cols=kwargs['ncols'])

        # for list of dataframes with one column (different index values)
        if style == 'multi-index':
            for series in data:
                ax.plot(series, label=series.columns[0], linewidth=0.25)
                ax.set_xlabel(data[0].index.name)
                ax.set_ylabel(kwargs['y_label'])
            Ploter.set_legend(ax, 4, font, cols=kwargs['ncols'])

        # annotations TODO: for multi styles
        if 'annot_texts' and'annot_points' in kwargs:
            for annotation in zip(kwargs['annot_points'], kwargs['annot_texts']):
                Ploter.set_annotation(ax, annotation[1], (annotation[0], data.loc[float(annotation[0]), data.columns.values[0]]))

        Ploter.set_ticks(ax, ticks)
        Ploter.set_fonts(ax, font)
        Ploter.save_figure(fig, title)

    @staticmethod
    def plot_multiple(style, data, title, size, font=12, ticks=5, **kwargs):

        # kwargs clean up
        kwargs = dict((k, v) for k, v in kwargs.items() if kwargs[k] is not None)

        fig, axes = Ploter.prepare_multiple_cx_figure(title, size, len(data))

        for d, ax in zip(data, axes):
            ax.plot(d)

        for ax in axes:
            Ploter.set_ticks(ax, ticks)
            Ploter.set_fonts(ax, font)

        Ploter.save_figure(fig, title)




