from modeller import *
from modeller.scripts import complete_pdb
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-n", "--nameTemplate")
parser.add_argument("-t", "--template")
parser.add_argument("-s", "--selectedNames", nargs='+')
args = parser.parse_args()

log.verbose()    # request verbose output
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

models=[args.nameTemplate+selected+'.pdb' for selected in args.selectedNames]
models.append(args.template)
profiles=['model_'+selected+'.profile' for selected in args.selectedNames]
profiles.append('template.profile')

for model, profile in zip(models, profiles):
    
    # read model file
    mdl = complete_pdb(env, model)

    # Assess with DOPE:
    s = selection(mdl)   # all atom selection
    s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=profile, normalize_profile=True, smoothing_window=15)
