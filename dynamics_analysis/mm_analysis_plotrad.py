# python 3
# script for plotting vmd radial distribution
# michaladammichalowski@gmail.com
# 25.11.15 - creation

import argparse
import numpy as np
import mm_lib_plots as mmplt


def read_rdf(rad_files):
    """
    :param rad_files: rdf files to read
    :return: list of rdf profiles
    """
    dist_rdf = []  # container for distances x
    dens_rdf = []  # container for densities y

    for rad_file in rad_files:
        rdf = open(rad_file).readlines()
        par_rdf = [[float(col) for col in row.split()] for row in rdf[0:-1]]  # first: create list rows
        dist_rdf.append([row[0] for row in par_rdf])  # second: select x columns from rows and add to container
        dens_rdf.append([row[1] for row in par_rdf])  # second: select y columns from rows and add to container

    return [dist_rdf, dens_rdf]


def plot_rdf(parseds_rdf):
    """
    :param parsed_xst: list of parsed rdf profiles
    :return: plots parsed rdf profiles (no return)
    """
    mmplt.plot_simple_multiple(parseds_rdf[0], parseds_rdf[1], "distance [A]", "density [rel]", args.labels, args.out_file)


parser = argparse.ArgumentParser()
parser.add_argument("-r", "--rad_files", nargs='+', help="rdf files to plot")
parser.add_argument("-l", "--labels", nargs='+', help="labels for data series")
parser.add_argument("-o", "--out_file", help="outfile name")
args = parser.parse_args()

rdfs = read_rdf(args.rad_files)
plot_rdf(rdfs)
