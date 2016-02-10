import subprocess as sp
import copy
import matplotlib.pyplot as plt
import argparse
import mm_lib_plots as mmplt

orange = '#f69855'
blue = '#3a81ba'
grey = '#666666'


def plot_profile(parseds_profs,  outfile):
    """
    :param parseds_profs: list of parsed rdf profiles
    :return: plots parsed rdf profiles (no return)
    """
    mmplt.plot_simple_multiple([parseds_profs[0], parseds_profs[0]], parseds_profs[1:3], \
                               "model []", "quality relative to best [rel]", ['molpdf', 'dope'], \
                               outfile, marker='.', linestyle='None', ylimit=[0, 1.1])

def md_table(title, data):
    """
    simple function to draw table summary
    :param title: name of the table
    :param data: list of model entries [model_number, molpdf, DOPE]
    :return: nothing
    """

    print('------------------------------------------------')
    print(title)
    print('------------------------------------------------')
    print('model    molpdf      DOPE')   
    for line in data:
        print('{0:5d}\t{1:7.0f}\t{2:9.0f}'.format(*line))
    print('------------------------------------------------')

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--logFile", dest="logFile", action="store")
parser.add_argument("-o", "--outFile", dest="outFile", action="store")
parser.add_argument("-n", "--modelNum", dest="modelNum", action="store", type=int)
args = parser.parse_args()


logFile = args.logFile
outFile = args.outFile
modelNum = args.modelNum


# read mark table from log into list of lists
modelList = [model.split() for model in sp.check_output(['grep', 'Summary of successfully produced models:', logFile, '-A', str(modelNum+2)]).decode('utf-8').split('\n')[3:-1]]

print(modelList)

# format each list of log
modelList = [[int(model[0][17:19]), float(model[1]), float(model[2])] for model in modelList]

# find best models
bestMolPDF = copy.deepcopy(min(modelList, key=lambda x: x[1]))
bestDOPE = copy.deepcopy(min(modelList, key=lambda x: x[2]))
print('Best molpdf: {}'.format(bestMolPDF))
print('Best DOPE: {}'.format(bestDOPE))

listByMOLPDF = sorted(modelList, key=lambda x: x[1])
listByDOPE = sorted(modelList, key=lambda x: x[2])

md_table('Best 25% by molpdf', listByMOLPDF[0:15])
md_table('Best 25% by DOPE', listByDOPE[0:15])

listBest = []
for model in listByMOLPDF[0:15]:
    if model in listByDOPE[0:15]:
        listBest.append(model)
listBestbyMOLPDF = sorted(listBest, key=lambda x: x[1])
md_table('Best common by molpdf', listBestbyMOLPDF)

# normalization (closer to one - better)
for model in modelList:
    model[1] = bestMolPDF[1]/model[1]
    model[2] = model[2]/bestDOPE[2]

# plots all models normalized (closer to one - better)
toPlot = list(zip(*modelList))
plot_profile(toPlot, args.outFile)
