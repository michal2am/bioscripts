# python 3
# parse modeller outifile
# michaladammichalowski@gmail.com
# ? - creation
# 18.02.16 - refactor
# EXAMPLE CALL: python3 mm_homology_eval_global.py -l mm_homology_create_gabar.log -o 3jae_global -n 256 -t 0.08

import subprocess as sp
import heapq
import copy
import argparse
import mm_lib_plots as mmplt


class ModelList:

    def __init__(self, log_file, model_num, treshold, out_file):
        """
        :param log_file: modeller log file
        :param model_num: number of models
        :param treshold: percent of top models to select
        :param out_file: file for plotting
        :return: none
        """

        self.log_file = log_file
        self.model_num = model_num
        self.tresh = treshold
        self.out_file = out_file

        self.best_molpdf = self.parse()[0]
        self.best_dope = self.parse()[1]
        self.best_common = self.parse()[2]
        self.normals = self.parse()[3]

    @staticmethod
    def get_extr(values, index):
        """
        :param values: line-like values iterable
        :param index: property to compare
        :return: lowest and highest line-like value
        """

        return heapq.nsmallest(1, values, key=lambda x: x[index])[0], \
               heapq.nlargest(1, values, key=lambda x: x[index])[0]

    @staticmethod
    def get_sort(values, index):
        """
        :param values: line-like values iterable
        :param index: property to compare
        :return: sorted line-like values iterable with position at [0]
        """

        return [[i+1] + model[:] for i, model in enumerate(sorted(values, key=lambda x: x[index]))]

    @staticmethod
    def get_top(values, perc):
        """
        :param values:
        :param perc:
        :return:
        """

        return values[0:int(len(values) * perc)]


    @staticmethod
    def normalized(values, index, bounds):
        """
        for properties 'the lower the better', 0-1 range 1 for best
        :param values: line-like values iterable
        :param index: properties to normalize [[i], ...]
        :param bounds: line-like values for boundaries [[best, worst], ...]
        :return: line-like values with normalized selected properties
        """

        for line in values:
            for ind in index:
                line[ind] = 1 - (line[ind] - bounds[ind-1][0][ind])/(bounds[ind-1][1][ind] - bounds[ind-1][0][ind])
        return values

    def parse(self):
        """
        :return: best by molpdf, best by dope, best common, all normalized
        """

        modelList = [model.split() for model in sp.check_output(['grep', 'Summary of successfully produced models:', \
                                                                 self.log_file, '-A', str(self.model_num+2)]) \
                                                                 .decode('utf-8').split('\n')[3:-1]]
        modelList = [[int(model[0][16:19]), float(model[1]), float(model[2])] for model in modelList]

        listByMOLPDF = self.get_sort(modelList, 1)
        listByDOPE = self.get_sort(modelList, 2)

        listBest = []
        for modelMol in self.get_top(listByMOLPDF, self.tresh):
            for modelDOPE in self.get_top(listByDOPE, self.tresh):
                if modelMol[1] == modelDOPE[1]:
                    pos = '{}/{}'.format(modelMol[0], modelDOPE[0])
                    listBest.append([pos] + modelMol[1:])

        bestMolPDF, worstMolPDF = self.get_extr(modelList, 1)
        bestDOPE, worstDOPE = self.get_extr(modelList, 2)

        tonorm = copy.deepcopy(modelList)
        listNorm = self.normalized(tonorm, [1, 2], [[bestMolPDF, worstMolPDF], [bestDOPE, worstDOPE]])

        return [listByMOLPDF, listByDOPE, listBest, listNorm]

    def print_table(self):
        """
        :return: prints table of best molpdf, dope and common
        """

        cols = '        model    molpdf      DOPE'
        sep = '-' * len(cols)
        for table, models in zip(['best by MOLPDF', 'best by DOPE', 'best common by MOLPDF'], \
                                 [self.get_top(self.best_molpdf, self.tresh), \
                                  self.get_top(self.best_dope, self.tresh), self.best_common]):
            print('{}\n{}\n{}\n{}\n{}'.format(sep, table, sep, cols, sep))
            for line in models:
                print('{0:>s}\t{1:5d}\t{2:7.0f}\t{3:9.0f}'.format(str(line[0]), *line[1:]))
            print(sep)

    def plot_profile(self):
        """
        :return: plots all models normalized scores
        """

        toplot = list(zip(*self.normals))
        mmplt.plot_simple_multiple([toplot[0], toplot[0]], toplot[1:3], \
                                   "model []", "quality relative to best [norm]", ['molpdf', 'dope'], \
                                   self.out_file, marker='.', linestyle='None', ylimit=[0, 1.1])


parser = argparse.ArgumentParser()
parser.add_argument("-l", "--logFile", dest="logFile", action="store")
parser.add_argument("-o", "--outFile", dest="outFile", action="store")
parser.add_argument("-n", "--modelNum", dest="modelNum", action="store", type=int)
parser.add_argument("-t", "--treshold", dest="treshold", action="store", type=float)
args = parser.parse_args()

models = ModelList(args.logFile, args.modelNum, args.treshold, args.outFile)
models.print_table()
models.plot_profile()
