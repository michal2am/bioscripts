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
    parseds_rdf = []

    for rad_file in rad_files:
        rdf = open(rad_file).readlines()
        par_rdf = [[float(col) for col in row.split()] for row in rdf[0:-1]]  # first: create list rows
        tog_rdf = [[row[0] for row in par_rdf], [row[1] for row in par_rdf]]  # second: select from rows
        # tog_rdf = [[float(row.split()[0]) for row in rdf[0:-1]], [float(row.split()[1]) for row in rdf[0:-1]]] # or select from string rows
        parseds_rdf.append(tog_rdf)

    return parseds_rdf


def plot_rdf(parseds_rdf):
    """
    :param parsed_xst: list of parsed rdf profiles
    :return: plots parsed rdf profiles (no return)
    """
    for parsed_rdf in parseds_rdf:
        mmplt.plot_simple(parsed_rdf[0], parsed_rdf[1], "distance [A]", "density [rel]", "equilibration_rdf")


parser = argparse.ArgumentParser()
parser.add_argument("-r", "--rad_files", nargs='+', help="rdf files to plot")
args = parser.parse_args()

rdfs = read_rdf(args.rad_files)
plot_rdf(rdfs)
