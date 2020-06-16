#import pylab
import modeller
import matplotlib.pyplot as plt

plot_range = 'all'

def md_plot(template, models_mod, xlab, ylab, models, plot_range, fontsize=10):

    fig, ax = plt.subplots()
    ax.set_color_cycle(['#4D4D4D','#5DA5DA', '#60BD68', '#F17CB0', '#B2912F', '#B276B2', '#DECF3F', '#F15854'])
    ax.tick_params(axis='both', which='major', labelsize=fontsize)

    if plot_range == 'all':
        ax.plot_euclidean(template, lw=0.6, color='#f69855', label='template')
        for model, number in zip(models_mod, models): 
            ax.plot_euclidean(model, lw=0.5, label='open hybrid model {0:02d}'.format(number))
    
    if plot_range == '1':
       ax.plot_euclidean(range(12, 350), template[674:1012], ls='-', color='#f69855', label='template')
       for model, number in zip(models_mod, models):
            ax.plot_euclidean(range(12, 350), model[0:338], lw=1, label='open hybrid model {0:02d}'.format(number))

    if plot_range == '2':
        ax.plot_euclidean(range(12, 350), template[0:338], ls='-', color='#f69855', label='template')
        for model, number in zip(models_mod, models):
            ax.plot_euclidean(range(12, 350), model[674:1012], lw=1, label='open hybrid model {0:02d}'.format(number))

    ax.locator_params(nbins=10)
    ax.set_xlabel(xlab, fontsize=fontsize)
    ax.set_ylabel(ylab, fontsize=fontsize)
    
    if plot_range == 'all': ax.set_xlim(0, 1700)
    else: ax.set_xlim(0, 370)
    ax.set_ylim(-0.08, -0.015)

    handles, labels = ax.get_legend_handles_labels()
   #lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=10)
    lgd = ax.legend(handles, labels, loc='lower right', ncol=2, fontsize=8)
    ax.grid('on')

    fig.set_size_inches(6, 4)
    fig.savefig('dope_profile_{}.png'.format(plot_range), dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight' )


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

e = modeller.environ()
a = modeller.alignment(e, file='gabar_MM.01chimeric.pir')

template = get_profile('chimeric.profile', a['chimeric'])

models = [1, 34, 53, 60, 42, 24, 49]
models_mod = []

for model in models:
    model_mod = get_profile('model_{0:02d}.profile'.format(model), a['gabar_MM.01'])
    models_mod.append(model_mod)

    profile_mod = open('dope_profile_model_{0:02d}.dat'.format(model), 'w')
    for i, res in enumerate(model_mod):
        profile_mod.write(str(i) + '   ' + str(res) + '\n')
    profile_mod.close()

# Plot the template and model profiles in the same plot for comparison:
#pylab.figure(1, figsize=(12,8))
#pylab.xlabel('Alignment position')
#pylab.ylabel('DOPE per-residue score')
#for model, number in zip(models_mod, models):
#   pylab.plot(model, linewidth=1, marker='o', ls='', label='open hybrid model {0:02d}'.format(number))
#pylab.plot(template, color='orange', linewidth=2, label='open hybrid template')
#pylab.legend()
#pylab.savefig('dope_profile.png', dpi=600)

print len(models_mod[0])
if plot_range == 'all': x_title = 'alignment position'
if plot_range == '1': x_title = 'residue number - first alpha subunit'
if plot_range == '2': x_title = 'residue number - second alpha subunit'
md_plot(template, models_mod, x_title, 'DOPE per residue score', models, plot_range)
#md_plot(template, models_mod, 'alignment position', 'DOPE per residue score', models)

