# python 3
# script for parsing namd outfiles
# michaladammichalowski@gmail.com
# 05.01.16 - creation
# EXAMPLE CALL:

import argparse
import numpy as np
import mm_lib_plots as mmplt


def read_out(template, selected_sims):
    """
    :param template: constant xst file name prefix
    :param selected_sims: iterative xst file name suffixes
    :return: concatenated xst file
    """
    out = open(template+selected_sims[0]+".out", 'r').readlines()

    for sim in selected_sims[1:]:
        out.extend(open(template+sim+".out", 'r').readlines()[2:])

    with open(template+"all.out", mode='wt', encoding='utf-8') as all_out:
        all_out.write(''.join(out))

    return out


def parse_out(out, timestep, save_freq):
    """
    :param out: concatenated xst file
    :param leaf_size: number of lipids in one leaf
    :param timestep: time step in fs
    :param save_freq: xst file save frequency
    :return: list of parsed xst parameters [time, apl, symmetry, total thickness, step number]
    """
    # time = 0
    # dtime = timestep*save_freq/1000000.

    key_out = ['TS', 'BOND', 'ANGLE', 'DIHED', 'IMPRP', 'ELECT', 'VDW', 'BOUNDARY', 'MISC', 'KINETIC', 'TOTAL', 'TEMP', \
               'POTENTIAL', 'TOTAL3', 'TEMPAVG', 'PRESSURE', 'GPRESSURE', 'VOLUME', 'PRESSAVG', 'GPRESSAVG']
    par_out = [{key: float(col) for (key, col) in zip(key_out, row.split()[1:])} \
               for row in out if len(row) > 1 and row.split()[0] == 'ENERGY:']  # first: create list rows

    # ts_all = [row[0]*timestep/10e5 for row in par_out]
    ts_all = [row['TS'] for row in par_out]
    return par_out


def plot_xst(par_out):
    """
    :param parsed_xst: list of parsed xst parameters
    :return: plots parsed xst parameters (no return)
    """
    mmplt.plot_simple([row['TS'] for row in par_out][500:], [row['POTENTIAL'] for row in par_out][500:], "step []", "potential energy [kcal/mol]", "potential", args.out_file + "_potential")
    mmplt.plot_simple([row['TS'] for row in par_out][500:], [row['TOTAL'] for row in par_out][500:], "step []", "total energy [kcal/mol]", "total", args.out_file + "_total")
    mmplt.plot_simple([row['TS'] for row in par_out][500:], [row['TEMP'] for row in par_out][500:], "step []", "temperature [K]", "temperature", args.out_file + "_temp")
    mmplt.plot_simple([row['TS'] for row in par_out][500:], [row['PRESSURE'] for row in par_out][500:], "step []", "pressure [bar]", "pressure", args.out_file + "_press")


parser = argparse.ArgumentParser()
parser.add_argument("-t", "--template", help="constant xst file name prefix")
parser.add_argument("-s", "--selected_sims", nargs='+', help="iterative xst file name suffixes")
parser.add_argument("-ts", "--timestep", type=int, help="time step in fs")
parser.add_argument("-sf", "--save_freq", type=int, help="xst file save frequency")
parser.add_argument("-o", "--out_file", help="outfile name")

args = parser.parse_args()

out = read_out(args.template, args.selected_sims)
par_out = parse_out(out, args.timestep, args.save_freq)
plot_xst(par_out)
