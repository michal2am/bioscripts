from modeller import *
import argparse
import pandas as pd
import mm_pandas_plot as mmpdplt

def r_enumerate(seq):
    """Enumerate a sequence in reverse order"""
    # Note that we don't use reversed() since Python 2.3 doesn't have it
    num = len(seq) - 1
    while num >= 0:
        yield num, seq[num]
        num -= 1

def get_profile(profile_file, seq):
    """Read `profile_file` into a Python array, and add gaps corresponding to
       the alignment sequence `seq`."""
    # Read all non-comment and non-blank lines from the file:
    f = open(profile_file).readlines()
    vals = []
    for line in f:
        if not line.startswith('#') and len(line) > 10:
            spl = line.split()
            vals.append(float(spl[-1]))
    # Insert gaps into the profile corresponding to those in seq:
    for n, res in r_enumerate(seq.residues):
        for gap in range(res.get_leading_gaps()):
            vals.insert(n, None)
    # Add a gap at position '0', so that we effectively count from 1:
    vals.insert(0, None)
    return vals


parser = argparse.ArgumentParser()
parser.add_argument("-a", "--aligFile")
parser.add_argument("-t", "--aligTemplate")
parser.add_argument("-s", "--aligSequence")
parser.add_argument("-m", "--selectedNames", nargs='+', type=int)
args = parser.parse_args()

e = environ()
a = alignment(e, file=args.aligFile)

template = get_profile('template.profile', a[args.aligTemplate])

models = args.selectedNames
models_mod = {}

for model in models:
    model_mod = get_profile('model_{0:02d}.profile'.format(model), a[args.aligSequence])
    model_mod.append(None)
    models_mod[model] = model_mod

    profile_mod = open('dope_profile_model_{0:02d}.dat'.format(model), 'w')
    for i, res in enumerate(model_mod):
        profile_mod.write(str(i) + '   ' + str(res) + '\n')
    profile_mod.close()

#print(template)
#model.append(None)
#print(model)
models_mod['template'] = template
all = pd.DataFrame(models_mod)
print(all)
for model in all.columns:
    all[model] = all[model] - all['template']


print(all)
ploter = mmpdplt.Ploter()
ploter.plot_single('single-index', all, 'fig_localleval', [10, 6], y_label='normalized score',
                  # axes_style={'ytvals': [-1, -0.2, -0.5, 0, 0.2, 0.5, 1], 'xtvals': [1000]},
                  lines_style={'linestyle': '-'},#, 'color': ['black', 'grey']}
                  legend_style={'loc': 'best', 'ncol': 5}
                   )

#md_plot(template, models_mod, 'alignment posiition', 'DOPE per residue score', models)