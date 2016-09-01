import MDAnalysis
import MDAnalysis.analysis.rms
import sys
import matplotlib.pyplot as plt

def md_plot (x, y, xlab, ylab, labe, savefile, fontsize=12):
    
    fig, ax = plt.subplots()
    orange = '#f69855'
    blue = '#3a81ba'
    grey = '#666666'

    ax.set_color_cycle([orange, blue, grey])

    ax.plot_euclidean(x, y[0], label=labe[0])
    ax.plot_euclidean(x, y[1], label=labe[1])
    ax.locator_params(nbins=3)
    ax.set_xlim([0, 120])
    ax.set_xlabel(xlab, fontsize=fontsize)
    ax.set_ylabel(ylab, fontsize=fontsize)

    handles, labels = ax.get_legend_handles_labels()
    #lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=10)
    lgd = ax.legend(handles, labels, loc='best', fontsize=10)
    ax.grid('on')

    fig.set_size_inches(4, 4)
    fig.savefig('fig_new/'+savefile+'.png', dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight' )


filename1 = sys.argv[1]
filename2 = sys.argv[2]
savefile = sys.argv[3]

#filename1 = 'Arg141_Tyr277_II.dat'
#filename2 = 'Arg141_Asp139_II.dat'
#savefile = 'beta_B6B7_M2M3_II'

cont1 = open(filename1).readlines()
cont2 = open(filename2).readlines()

print cont2
labe = [cont1[0], cont2[0]]

nowtime = 0.0
dtime = 50*(1.0*25000/1000000) 
time  = []
dist1 = []
dist2 = []

for con1, con2 in zip(cont1[1:], cont2[1:]):
	con1 = con1.split()
	con2 = con2.split()
	time.append(nowtime)
	dist1.append(con1[1])
	dist2.append(con2[1])
	nowtime += dtime

md_plot(time, [dist1, dist2], 'time [ns]', 'distance [A]', labe, savefile)
