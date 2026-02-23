import subprocess as sp
import matplotlib.pyplot as plt
import numpy as np

orange = '#f69855'
blue = '#3a81ba'
grey = '#666666'


def md_plot(data, xlab, ylab, labe, fontsize=12):

    fig, ax = plt.subplots()
    ax.set_color_cycle(['#4D4D4D','#5DA5DA', '#60BD68', '#F17CB0', '#B2912F', '#B276B2', '#DECF3F', '#F15854'])
    pore_range_temp=[120, -550]
    pore_range=[145, -550]
    
    for (series, model) in zip(data, labe):
        series = list(zip(*series))
        if model == 'chimeric': 
           ax.plot_euclidean(series[1][pore_range_temp[0]:pore_range_temp[1]], series[0][pore_range_temp[0]:pore_range_temp[1]], label='template', color='#f69855')
        else: 
           lab = 'closed hybrid model {0:02d}'.format(model)
           ax.plot_euclidean(series[1][pore_range[0]:pore_range[1]], series[0][pore_range[0]:pore_range[1]], label=lab, lw=1)
 
    #ax.locator_params(nbins=20)
    ax.set_xlabel(xlab, fontsize=fontsize)
    ax.set_ylabel(ylab, fontsize=fontsize)
   #ax.set_xlim([0, 120])
   #ax.set_ylim([0.75, 1.02])

    handles, labels = ax.get_legend_handles_labels()
    #lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=10)
    lgd = ax.legend(handles, labels, loc='best', fontsize=10)
    ax.grid('on')

    fig.set_size_inches(7.5, 6)
    fig.savefig('all_profiles.png', dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight' )

def md_table(title, data):

    print('------------------------------------------------')
    print(title)
    print('------------------------------------------------')
    print('model    radius [A]')
    for model in data.keys():
        print('{0}\t{1:7.2f}'.format(model, data[model]))
    print('------------------------------------------------')


def position_model(models, path):

    for model in models:

        if model == 'chimeric':
            infile = 'chimeric_fit.pdb'
            outfile = 'chimeric_fit_pos.pdb'
        else:
            infile = 'gabar_MM.01.B999900{0:02d}_fit.pdb'.format(model)
            outfile = 'model_{0:02d}_fit_pos.pdb'.format(model)
   
        infile = path  + infile 
     
        vmd_script = 'mol load pdb {0} \nsource orienter.tcl\nsource centerer.tcl\n[atomselect top all] move [trans x 180]\n[atomselect top all] writepdb {1}\nmol delete top\nresetpsf'.format(infile, outfile)
        sp.check_output(['vmd', '-e'], input=vmd_script, universal_newlines=True)
    
def create_profile(models):

    for model in models:

        if model == 'chimeric':
            infile = 'chimeric_fit_pos.pdb'
            outfile = 'chimeric_fit_pos.profile'
        else:
            infile = 'model_{0:02d}_fit_pos.pdb'.format(model)
            outfile = 'model_{0:02d}_fit_pos.profile'.format(model)

        hole2_script = 'coord {0} \nradius /home/michal_am/Copy/MD/software_packages/hole2/rad/simple.rad \ncvect 0.0 0.0 1.0 \ncpoint 0.  0.  0.'.format(infile)
        profile = open(outfile, 'w')
        profile.write(sp.check_output(['/opt/hole2/exe/hole'], input=hole2_script, universal_newlines=True))
        profile.close()   

def read_profile(models):
    
    min_rad = {}
    profiles = []

    for model in models:

        if model == 'chimeric':
            infile = 'chimeric_fit_pos.profile'
            #outfile = 'chimeric_fit_pos.profileGrph'
        else:
            infile = 'model_{0:02d}_fit_pos.profile'.format(model)
            #outfile = 'model_{0:02d}_fit_pos.profileGrph'.format(model)

        prof = open(infile, 'r').readlines()
        profile = []
    
        for line in prof:
            line = line.split()
            if len(line) == 5 and line[4] in ['(sampled)', '(mid-point)']:
                profile.append([float(col) for col in line[0:2]])
            if len(line) == 5 and line[0] == 'Minimum':
                min_rad[model] =  float(line[3])
        profiles.append(profile) 
     
    md_table('Minimum pore radiuses', min_rad)

    md_plot(profiles, 'pore radius [A]', 'distance along pore axis [A]', models)

models = ['chimeric', 1, 34, 53, 60, 42, 24, 49]
path = ''

#position_model(models, path)
create_profile(models)
read_profile(models)

