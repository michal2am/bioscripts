# python 3
# plots histogram of vmd inline coordinates
# michaladammichalowski@gmail.com
# 11.04.16 - creation
# EXAMPLE CALL:

import argparse
import numpy as np
import mm_lib_plots as mmplt


class Coordinates:

    def __init__(self, infile, outfile, label):
        """
        :param infile:
        :param outfile:
        :return:
        """

        self.infile = infile
        self.outfile = outfile
        self.label = label

        self.coors = self.read_coors()

    def read_coors(self):
        """
        :param args:
        :return:
        """
        coordinates = []
        for coor_file in self.infile:
            with open(coor_file, 'r') as coors:
                coors = coors.readlines()
                coors = [[float(col) for col in row.split()] for row in coors]  # first: create list rows
                coors = [item for sublist in coors for item in sublist]       # second: put rows alltogether
                coors = [ coor for coor in coors if (-25 < coor < 35)]
                coordinates.append(coors)
        return coordinates

    def plot_coors(self):
        """
        :param parseds_thc:
        :return:
        """
        mmplt.plot_histogram(self.coors, "density [norm]", "pore axis coordinate [A]", self.label, self.outfile,
                             sizex=1.5, sizey=1.5, fontsize=6, ranges=True)


parser = argparse.ArgumentParser()
parser.add_argument("-c", "--infile", nargs='+', help="coordinate to plot")
parser.add_argument("-p", "--outfile", help="file to save plot")
parser.add_argument("-l", "--label", nargs='+', help="plot label")
args = parser.parse_args()

water_coors = Coordinates(args.infile, args.outfile, args.label)
water_coors.plot_coors()
