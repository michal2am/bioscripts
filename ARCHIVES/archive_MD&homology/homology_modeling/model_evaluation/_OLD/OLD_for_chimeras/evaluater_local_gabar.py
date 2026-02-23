from modeller import *
from modeller.scripts import complete_pdb

log.verbose()    # request verbose output
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

models = ['chimeric', 1, 34, 53, 60, 42, 24, 49]

for model in models:
    
    if model == 'chimeric':
        infile = 'chimeric_fit.pdb'
        outfile = 'chimeric.profile'
    else:
        infile = 'gabar_MM.01.B999900{0:02d}_fit.pdb'.format(model)
        outfile = 'model_{0:02d}.profile'.format(model)

    # read model file
    mdl = complete_pdb(env, infile)

    # Assess with DOPE:
    s = selection(mdl)   # all atom selection
    s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=outfile, normalize_profile=True, smoothing_window=15)
