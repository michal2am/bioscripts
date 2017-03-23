import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt


class Ploter:

    def __init__(self):

        # seaborn globals

        sns.set_context('paper')
        sns.set_style('ticks', {'legend.frameon': True})
        sns.set_palette("Set1", 8)

    @classmethod
    def set_rc(cls, fontsize):
        """
        sets font sizes via rcParam, use before figure generation
        :param fontsize: font size
        """

        mpl.rcParams['font.family'] = 'sans-serif'                              # global
        mpl.rcParams['font.size'] = fontsize                                    # only for relative settings

        mpl.rcParams['axes.labelsize']  = fontsize                              # element specific font size
        mpl.rcParams['axes.titlesize']  = fontsize
        mpl.rcParams['xtick.labelsize'] = fontsize
        mpl.rcParams['ytick.labelsize'] = fontsize
        mpl.rcParams['legend.fontsize'] = fontsize

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
        if 'ytvals' in kwargs: ax.yaxis.set_ticks(kwargs['ytvals'])

        if 'xlim' in kwargs: ax.set_xlim(kwargs['xlim'])
        if 'ylim' in kwargs: ax.set_ylim(kwargs['ylim'])

    @classmethod
    def set_legend(cls, ax, **kwargs):
        """
        sets legend parameters
        :param ax: axes
        """

        leg = ax.legend(**kwargs)
        for legobj in leg.legendHandles:
            legobj.set_linewidth(2.0)

    @staticmethod
    def plot_single(style, data, title, size, kind='line', font=20, **kwargs):
        """
        plot single plot of one dataseries
        :param style: 'single' or 'multi-indepx'
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
        if style == 'single-index':

            data.plot(ax=ax, kind=kind, **kwargs.get('lines_style', {}))
            ax.set_xlabel(data.index.name)
            if len(data.columns.values) == 1:
                ax.set_ylabel(data.columns.values[0])
            else:
                ax.set_ylabel(kwargs['y_label'])
                Ploter.set_legend(ax, **kwargs.get('legend_style', {}))

        Ploter.set_axes(ax, **kwargs.get('axes_style', {}))
        Ploter.save_figure(fig, title, kwargs.get('rect', (0, 0, 1, 1)))
