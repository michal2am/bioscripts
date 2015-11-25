#import pylab
import modeller
import matplotlib.pyplot as plt
import argparse


def md_plot(template, models_mod, xlab, ylab, models, fontsize=10):

    fig, ax = plt.subplots()
    ax.set_color_cycle(['#4D4D4D','#5DA5DA', '#60BD68', '#F17CB0', '#B2912F', '#B276B2', '#DECF3F', '#F15854'])
    ax.tick_params(axis='both', which='major', labelsize=fontsize)

    ax.plot(template, lw=0.6, color='#f69855', label='template')
    for model, number in zip(models_mod, models): 
        ax.plot(model, lw=0.5, label='open hybrid model {0:02d}'.format(number))
    
    ax.locator_params(nbins=10)
    ax.set_xlabel(xlab, fontsize=fontsize)
    ax.set_ylabel(ylab, fontsize=fontsize)
    
    # ax.set_xlim(0, 1700)
    # ax.set_ylim(-0.08, -0.015)

    handles, labels = ax.get_legend_handles_labels()
   #lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=10)
    lgd = ax.legend(handles, labels, loc='lower right', ncol=2, fontsize=8)
    ax.grid('on')

    fig.set_size_inches(6, 4)
    fig.savefig('dope_profile.png', dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight' )


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
    f = file(profile_file)
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

e = modeller.environ()
a = modeller.alignment(e, file=args.aligFile)

template = get_profile('template.profile', a[args.aligTemplate])

models = args.selectedNames
models_mod = []

for model in models:
    model_mod = get_profile('model_{0:02d}.profile'.format(model), a[args.aligSequence])
    models_mod.append(model_mod)

    profile_mod = open('dope_profile_model_{0:02d}.dat'.format(model), 'w')
    for i, res in enumerate(model_mod):
        profile_mod.write(str(i) + '   ' + str(res) + '\n')
    profile_mod.close()

md_plot(template, models_mod, 'alignment posiition', 'DOPE per residue score', models)

