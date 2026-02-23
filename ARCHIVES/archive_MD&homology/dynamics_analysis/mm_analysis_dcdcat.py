# python 2
# script for joining dcd files
# michaladammichalowski@gmail.com
# xx.xx.14 - creation
# 26.11.15 - refactor

from MDAnalysis import *
import argparse


def cat_dcd(template, selected_sims, timestep, save_freq, psf_file):
        """
        :param template: constant dcd file name prefix
        :param selected_sims: iterative dcd file name suffixes
        :param timestep: time step in fs
        :param save_freq: xst file save frequency
        :param psf_file: full psf file name
        :return: writes all trajectories to one dcd (nothing)
        """
        dcd = list()
        fs2akma = 0.0204582651
        efts = timestep*fs2akma

        for traj in selected_sims:
            dcd.append(template+traj+".dcd")

        print "Merging: " + ', '.join(dcd)
        u = Universe(psf_file, dcd)
        w = Writer(template+"_all.dcd", u.atoms.n_atoms, delta=efts, step=save_freq )

        for ts in u.trajectory:
            w.write(u)

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--template", help="constant dcd file name prefix")
parser.add_argument("-s", "--selected_sims", nargs='+', help="iterative dcd file name suffixes")
parser.add_argument("-ts", "--timestep", type=int, help="time step in fs")
parser.add_argument("-sf", "--save_freq", type=int, help="xst file save frequency")
parser.add_argument("-psf", "--psf_file", help="full psf file name")
args = parser.parse_args()

cat_dcd(args.template, args.selected_sims, args.timestep, args.save_freq, args.psf_file)