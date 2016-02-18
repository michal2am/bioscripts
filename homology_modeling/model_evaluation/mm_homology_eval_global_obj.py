# python 3
# parse modeller outifile
# michaladammichalowski@gmail.com
# ? - creation
# 18.02.16 - refactor
# EXAMPLE CALL:

import subprocess as sp
import copy
import argparse
import mm_lib_plots as mmplt


class ModelList:

    def __init__(self, log_file, model_num, out_file):

        self.log_file = log_file
        self.model_num = model_num
        self.out_file = out_file

        self.best_molpdf = self.get_best()[0]
        self.best_dope = self.get_best()[1]
        # self.best_common = self.get_best()[2]

    @staticmethod
    def get_extr(values, index):
        return [copy.deepcopy(min(values, key=lambda x: x[index])), copy.deepcopy(max(values, key=lambda x: x[index]))]

    @staticmethod
    def get_sort(values, index):
        sort = sorted(values, key=lambda x: x[index])
        return [[i, model[:]] for i, model in enumerate(sort)]

    def get_best(self):

        modelList = [model.split() for model in sp.check_output(['grep', 'Summary of successfully produced models:', self.log_file, '-A', str(self.model_num+2)]).decode('utf-8').split('\n')[3:-1]]
        modelList = [[int(model[0][16:19]), float(model[1]), float(model[2])] for model in modelList]

        bestMolPDF, worstMolPDF = self.get_extr(modelList, 1)
        bestDOPE, worstDOPE = self.get_extr(modelList, 2)

        listByMOLPDF = self.get_sort(modelList, 1)
        listByDOPE = self.get_sort(modelList, 2)

        listBest = []
        for model in listByMOLPDF[0:15]:
            if model in listByDOPE[0:15]:
                listBest.append(model)
        listBestbyMOLPDF = sorted(listBest, 1)

        return [listByMOLPDF, listByDOPE]

    def normalize(self):



    def print_table(self, title):

        print('------------------------------------------------')
        print(title)
        print('------------------------------------------------')
        print('no     model    molpdf      DOPE')
        for line in []:
            print('{0:5d}\t{1:7.0f}\t{2:9.0f}'.format(*line))
            print('------------------------------------------------')


parser = argparse.ArgumentParser()
parser.add_argument("-l", "--logFile", dest="logFile", action="store")
parser.add_argument("-o", "--outFile", dest="outFile", action="store")
parser.add_argument("-n", "--modelNum", dest="modelNum", action="store", type=int)
args = parser.parse_args()

models = ModelList(args.logFile, args.modelNum, args.outFile)
print(models.best_molpdf[0:5])