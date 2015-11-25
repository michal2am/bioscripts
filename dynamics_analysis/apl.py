# script for calculating area per lipid and system thickness
# michaladammichalowski@gmail.com
# xx.xx.14 - creation
# 25.11.15 - refactor

import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
# mpl.use('Agg')

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--template", help="constant xst file name prefix")
parser.add_argument("-s", "--selected_sims", nargs='+', help="iterative xst file name suffixes")
parser.add_argument("-l", "--leaf_size", type=int, help="number of lipids in one leaf")
parser.add_argument("-ts", "--timestep", type=int, help="time step in fs")
parser.add_argument("-sf", "--save_freq", type=int, help="xst file save frequency")
args = parser.parse_args()


def read_xst(template, selected_sims):
    """
    :param template: constant xst file name prefix
    :param selected_sims: iterative xst file name suffixes
    :return: concatenated xst file
    """
    xst = open(template+selected_sims[0]+".xst", 'r').readlines()
    for sim in selected_sims[1:]:
        xst.extend(open(template+sim+".xst", 'r').readlines()[2:])

    with open(template+"all.xst", mode='wt', encoding='utf-8') as all_xst:
        all_xst.write(''.join(xst))
    return xst

xst = read_xst(args.template, args.selected_sims)



def parse_xst(xst, leaf_size, timestep, save_freq):
    """
    :param xst: concatenated xst file
    :param leaf_size: number of lipids in one leaf
    :param timestep: time step in fs
    :param save_freq: xst file save frequency
    :return:
    """
    # time = 0
    # dtime = timestep*save_freq/1000000.

    par_xst = [[float(col) for col in row.split()] for row in xst[2:]]

    timeAll = [row[0]*timestep*10e-6 for row in par_xst]
    apl_all = [row[1]*row[5]/leaf_size for row in par_xst]
    sym_all = [row[1]/row[5] for row in par_xst]
    thc_all = [row[9] for row in par_xst]

parse_xst(xst, args.leaf_size, args.timestep, args.save_freq)

"""
fig = plt.figure(figsize=(6, 8))
ax = fig.add_subplot(311)
ax.plot(timeAll, aplAll, 'k-',  label="area per lipid")
ax.legend(loc="best")
ax.set_ylabel("APL [A]")
ax.set_xlabel("time [ns]")
ay = fig.add_subplot(312)
ay.plot(timeAll, thcAll, 'k--', label="system thickness")
ay.legend(loc="best")
ay.set_ylabel("thickness [A]")
ay.set_xlabel("time [ns]")
az = fig.add_subplot(313)
az.plot(timeAll, symAll, 'k--', label="system symmetry")
az.legend(loc="best")
az.set_ylabel("symmetry X:Y")
az.set_xlabel("time (ps)")

fig.savefig("apl.pdf")
"""