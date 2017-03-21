# python 3
# parse modeller outifile
# michaladammichalowski@gmail.com
# ? - creation
# 18.02.16 - refactor
# 21.03.17 - refactor
# EXAMPLE CALL: python3 mm_homology_eval_global.py -l mm_homology_create_gabar.log -o 3jae_global -n 256 -t 0.08

import subprocess as sp
import argparse
import numpy as np
import pandas as pd
import mm_pandas_plot as mm_plt

class ModelList:

    def __init__(self, log_file, model_num, treshold, out_file):
        """
        :param log_file: modeller log file
        :param model_num: number of models
        :param treshold: percent of top models to select
        :param out_file: file for plotting
        """

        self.log_file, self.model_num, self.tresh, self.out_file = log_file, model_num, int(treshold*model_num), out_file
        self.best_molpdf, self.best_dope, self.best_common, self.normals = self.parse()

    def parse(self):
        """
        :return: best by molpdf, best by dope, best common, best common normalized
        """

        modelList = [model.split() for model in sp.check_output(['grep', 'Summary of successfully produced models:', \
                                                                 self.log_file, '-A', str(self.model_num+2)]) \
                                                                 .decode('utf-8').split('\n')[3:-1]]
        modelList = [[int(model[0][15:19]), float(model[1]), float(model[2])] for model in modelList]

        by_model = pd.DataFrame(modelList, columns=['model', 'molpdf', 'dope']).set_index(keys='model', drop=True)

        by_molpdf = by_model.sort_values(by=['molpdf']).iloc[0:self.tresh]
        by_dope = by_model.sort_values(by=['dope']).iloc[0:self.tresh]

        common = by_molpdf.index.intersection(by_dope.index)
        c_by_molpdf = by_model.loc[common].sort_values(by=['molpdf'])
        n_by_model = by_model.apply(lambda x: (x - np.mean(x))/(np.min(x) - np.max(x)))

        print('Top by molpdf:')
        print(by_molpdf)
        print('Top by dope')
        print(by_dope)
        print('Top common:')
        print(c_by_molpdf)

        return by_molpdf, by_dope, c_by_molpdf, n_by_model

    def plot_profile(self):
        """
        :return: plots all models normalized scores
        """
        ploter = mm_plt.Ploter()
        ploter.plot_single('single-index', self.normals, 'fig_globaleval', [6, 6], y_label='normalized score',
                           axes_style={'ytvals': [-1, -0.2, -0.5, 0, 0.2, 0.5, 1], 'xtvals': [1000]},
                           lines_style={'linestyle': 'None', 'marker': 'o', 'color': ['black', 'grey']},
                           legend_style={'loc': 'best', 'ncol': 2}
                           )

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--logFile", dest="logFile", action="store")
parser.add_argument("-o", "--outFile", dest="outFile", action="store")
parser.add_argument("-n", "--modelNum", dest="modelNum", action="store", type=int)
parser.add_argument("-t", "--treshold", dest="treshold", action="store", type=float)
args = parser.parse_args()

models = ModelList(args.logFile, args.modelNum, args.treshold, args.outFile)
models.plot_profile()
