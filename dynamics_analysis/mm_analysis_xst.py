# python 3
# script for calculating area per lipid and system thickness
# michaladammichalowski@gmail.com
# xx.xx.14 - creation
# 25.11.15 - refactor

import argparse
import numpy as np
import mm_lib_plots as mmplt


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


def parse_xst(xst, leaf_size, timestep, save_freq):
    """
    :param xst: concatenated xst file
    :param leaf_size: number of lipids in one leaf
    :param timestep: time step in fs
    :param save_freq: xst file save frequency
    :return: list of parsed xst parameters [time, apl, symmetry, total thickness, step number]
    """
    # time = 0
    # dtime = timestep*save_freq/1000000.

    par_xst = [[float(col) for col in row.split()] for row in xst[2:]]

    tim_all = [row[0]*timestep*10e-6 for row in par_xst]
    ste_all = [row[0] for row in par_xst]
    apl_all = [row[1]*row[5]/leaf_size for row in par_xst]
    sym_all = [row[1]/row[5] for row in par_xst]
    thc_all = [row[9] for row in par_xst]

    print("Area per lipid: AVR: {0}, SD: {1}".format(np.mean(apl_all), np.std(apl_all)))
    print("Total thickness: AVR: {0}, SD: {1}".format(np.mean(thc_all), np.std(thc_all)))
    print("Patch symmetry: AVR: {0}, SD: {1}".format(np.mean(sym_all), np.std(sym_all)))

    return [tim_all, apl_all, sym_all, thc_all, ste_all]


def plot_xst(parsed_xst):
    """
    :param parsed_xst: list of parsed xst parameters
    :return: plots parsed xst parameters (no return)
    """
    mmplt.plot_simple(parsed_xst[4], parsed_xst[1], "step", "APL [A]", "equilibration_apl")
    mmplt.plot_simple(parsed_xst[4], parsed_xst[2], "step", "symmetry", "equilibration_symmetry")
    mmplt.plot_simple(parsed_xst[4], parsed_xst[3], "step", "total thickness [A]", "equilibration_totthc")

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--template", help="constant xst file name prefix")
parser.add_argument("-s", "--selected_sims", nargs='+', help="iterative xst file name suffixes")
parser.add_argument("-l", "--leaf_size", type=int, help="number of lipids in one leaf")
parser.add_argument("-ts", "--timestep", type=int, help="time step in fs")
parser.add_argument("-sf", "--save_freq", type=int, help="xst file save frequency")
args = parser.parse_args()

xst = read_xst(args.template, args.selected_sims)
parsed_xst = parse_xst(xst, args.leaf_size, args.timestep, args.save_freq)
plot_xst(parsed_xst)