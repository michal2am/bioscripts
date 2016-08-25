import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import seaborn as sns
import pandas as pd


class ABF_results:

    def __init__(self, filename):
        """
        reads abf namd output files
        :param filename: output filename without extension
        :return: none
        """

        self.filename = filename
        self.fcount, self.fgrad, self.fpmf = [self.filename + ext for ext in ['.count', '.grad', '.pmf']]
        self.count, self.grad, self.pmf = self.read_result(self.fcount, 3, ['z [A]', 'count']),\
                                          self.read_result(self.fgrad, 3, ['z [A]', 'grad']),\
                                          self.read_result(self.fpmf, 1, ['z [A]', 'pmf'])

    @staticmethod
    def read_result(fname, skip, header):
        """
        read space separated data file
        :param fname: file name
        :param skip: number of top rows to skip
        :param header: header (list of strings) for dataframe
        :return: dataframe with parsed data
        """

        data = pd.read_csv(fname, delim_whitespace=True, skiprows=skip, names=header, index_col=0)
        return data

    @staticmethod
    def plot_vertical_result(data, title):
        """
        plots mutliple dataframes in vertical set
        :param data: dataframe top plot
        :param title: title of graph
        :return: none
        """

        dnum = len(data)

        sns.set_style("ticks", {'legend.frameon': True})
        sns.set_context("poster")
        sns.set_palette('Paired', dnum)

        fig = plt.figure()
        fig.canvas.set_window_title(title)
        gs = grd.GridSpec(dnum, 1, hspace=0.5)

        for i, d in enumerate(data):
            ax = plt.subplot(gs[i, 0])
            d.plot(ax=ax)
            fig.add_subplot(ax)

        sns.despine()
        plt.show()

filename = 'step6.6ddabf_equilibration'
title = 'ABF results'

results = ABF_results(filename)
results.plot_vertical_result([results.count, results.grad, results.pmf], title)

