# python 3
# parse modeller outifile
# michaladammichalowski@gmail.com
# ? - creation
# 18.02.16 - refactor
# EXAMPLE CALL:

import subprocess as sp
import heapq
import argparse
import mm_lib_plots as mmplt


class ModelList:

    def __init__(self, log_file, model_num, out_file):
        """
        :param log_file: modeller log file
        :param model_num: number of models
        :param out_file: file for plotting
        :return: none
        """

        self.log_file = log_file
        self.model_num = model_num
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

        modelList = [model.split() for model in sp.check_output(['grep', 'Summary of successfully produced models:', self.log_file, '-A', str(self.model_num+2)]).decode('utf-8').split('\n')[3:-1]]
        modelList = [[int(model[0][16:19]), float(model[1]), float(model[2])] for model in modelList]

        listByMOLPDF = self.get_sort(modelList, 1)
        listByDOPE = self.get_sort(modelList, 2)

        listBest = []
        for modelMol in listByMOLPDF[0:int(self.model_num * 0.05)]:
            for modelDOPE in listByDOPE[0:int(self.model_num * 0.05)]:
                if modelMol[1] == modelDOPE[1]:
                    pos = '{}/{}'.format(modelMol[0], modelDOPE[0])
                    listBest.append([pos] + modelMol[1:])

        bestMolPDF, worstMolPDF = self.get_extr(modelList, 1)
        bestDOPE, worstDOPE = self.get_extr(modelList, 2)
        listNorm = self.normalized(modelList, [1, 2], [[bestMolPDF, worstMolPDF], [bestDOPE, worstDOPE]])

        return [listByMOLPDF, listByDOPE, listBest, listNorm]

    def print_table(self, title):

        print('-'*40)
        print(title)
        print('-'*40)
        print('no     model    molpdf      DOPE')
        print('-'*40)
        for line in []:
            print('{0:5d}\t{1:7.0f}\t{2:9.0f}'.format(*line))
            print('-'*40)
'''
def plot_profile(parseds_profs,  outfile):
    """
    :param parseds_profs: list of parsed rdf profiles
    :return: plots parsed rdf profiles (no return)
    """
    mmplt.plot_simple_multiple([parseds_profs[0], parseds_profs[0]], parseds_profs[1:3],
                               "model []", "quality relative to best [norm]", ['molpdf', 'dope'], \
                               outfile, marker='.', linestyle='None', ylimit=[0, 1.1])
'''
parser = argparse.ArgumentParser()
parser.add_argument("-l", "--logFile", dest="logFile", action="store")
parser.add_argument("-o", "--outFile", dest="outFile", action="store")
parser.add_argument("-n", "--modelNum", dest="modelNum", action="store", type=int)
args = parser.parse_args()

models = ModelList(args.logFile, args.modelNum, args.outFile)
print(models.best_molpdf[0:5])
print(models.best_dope[0:5])
print(models.best_common[0:5])
print(models.normals[0:5])

'''
# plots all models normalized (closer to one - better)
toPlot = list(zip(*modelList))
plot_profile(toPlot, args.outFile)
'''
