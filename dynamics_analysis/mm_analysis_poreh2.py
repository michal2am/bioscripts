# python 3
# plots pore profiles using hole2
# michaladammichalowski@gmail.com
# 16.03.15 - creation
# EXAMPLE CALL:

import argparse as ap
import os
import subprocess as sp
import numpy as np
import mm_lib_plots as mmplt


# hole2 fixed file paths
hole2 = '/opt/hole2/exe/hole'
hole2_rads = '/opt/hole2/rad/simple.rad'


class PoreModel:

    def __init__(self, pdb_dir, name):
        self.pdb_dir = pdb_dir
        self.name = name

        self.positioned = self.get_positions()
        self.profiles = self.get_profile()

    def set_position(self):
        """
        calls complete_positioner to align and center models
        :return:
        """
        sp.check_call(['complete_positioner.sh', self.pdb_dir])

    def get_positions(self):
        """
        :return:
        """
        return [pdb for pdb in os.listdir(self.pdb_dir + '/positioned') if pdb.endswith('_pos.pdb')]

    def get_profile(self):
        """
        :return:
        """
        profiles = []

        for pdb in self.positioned:
            outfile = os.path.splitext(pdb)[0] + '.dat'
            hole2_script = 'coord {0} \nradius {1} \ncvect 0.0 0.0 1.0 \ncpoint 0.  0.  -15.'\
                .format(pdb, hole2_rads)
            with open(outfile, 'w') as pore_file:
                pore_file.write(sp.check_output([hole2], stdin=hole2_script, universal_newlines=True))
                profiles.append(pore_file)

        return profiles

    def parse_profile(self):
        """
        :return:
        """
        parsed_profiles = []

        for profile in self.profiles:
            with open(profile, 'r').readlines() as prof:
                profile = []
                for line in prof:
                    line = line.split()
                    if len(line) == 5 and line[4] in ['(sampled)', '(mid-point)']:
                        profile.append([float(col) for col in line[0:2]])
                    if len(line) == 5 and line[0] == 'Minimum':
                        min_rad[model] =  float(line[3])
                profiles.append(profile)




parser = ap.ArgumentParser()
parser.add_argument("-d", "--pdb_dirs", nargs='+', help="directories with pdb files to analyze")
parser.add_argument("-n", "--pdb_names", nargs='+', help="model names")
args = parser.parse_args()
