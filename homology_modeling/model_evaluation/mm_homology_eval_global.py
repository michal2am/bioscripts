


def plot_profile(parseds_profs,  outfile):
    """
    :param parseds_profs: list of parsed rdf profiles
    :return: plots parsed rdf profiles (no return)
    """
    mmplt.plot_simple_multiple([parseds_profs[0], parseds_profs[0]], parseds_profs[1:3],
                               "model []", "quality relative to best [norm]", ['molpdf', 'dope'], \
                               outfile, marker='.', linestyle='None', ylimit=[0, 1.1])

def md_table(title, data):
    """
    simple function to draw table summary
    :param title: name of the table
    :param data: list of model entries [model_number, molpdf, DOPE]
    :return: nothing
    """



md_table('Best 25% by molpdf', listByMOLPDF[0:15])
md_table('Best 25% by DOPE', listByDOPE[0:15])


md_table('Best common by molpdf', listBestbyMOLPDF)

# normalization (closer to one - better)
for model in modelList:
    model[1] = 1 - (model[1] - bestMolPDF[1])/(worstMolPDF[1] - bestMolPDF[1])
    model[2] = 1 - (model[2] - bestDOPE[2])/(worstDOPE[2] - bestDOPE[2])

# plots all models normalized (closer to one - better)
toPlot = list(zip(*modelList))
plot_profile(toPlot, args.outFile)
