# python 3
# plots pore profiles using hole2
# michaladammichalowski@gmail.com
# 16.03.15 - creation
# EXAMPLE CALL:

import argparse as ap
import os
import subprocess as sp
import logging as lg
import numpy as np
import mm_lib_plots as mmplt


# hole2 fixed file paths
hole2 = '/opt/hole2/exe/hole'
hole2_rads = '/opt/hole2/rad/simple.rad'
positioner = '/home/mm/Pycharm/bioscripts/homology_modeling/model_positioning/complete_positioner.sh'

class PoreModel:

    def __init__(self, pdb_dir, name):
        self.pdb_dir = pdb_dir
        self.name = name
        self.labels = self.set_names()

        # self.set_position()
        # self.get_profile()
        self.parsed_profiles, self.min_Rads = self.parse_profiles()

    @staticmethod
    def find_files(directory, extension):
        """
        :param directory: starting directory
        :param extension: files to look for
        :return: list of relative paths to files
        """
        file_list = []
        for dir_, _, files in os.walk(directory):
            for file in files:
                if file.endswith(extension):
                    relDir = os.path.relpath(dir_, directory)
                    file_list.append(os.path.join(relDir, file))
        return file_list

    def set_names(self):
        """
        :return:
        """
        #return [os.path.splitext(os.path.split(pdb)[1])[0] for pdb in self.find_files(self.pdb_dir, '.pdb')]
        return [str(i) for i in range(12)]

    def set_position(self):
        """
        calls complete_positioner to align and center models
        :return:
        """
        log.info('Searching in dirs: {} for models {}'.format(self.pdb_dir, self.name))
        log.info('Positioning models according to principal axes and geometry center')
        sp.check_call([positioner, self.pdb_dir])

    def get_profile(self):
        """
        :return:
        """
        log.info('Reading positioned models from {}'.format(self.pdb_dir))
        log.info('Esitmating pore profiles with Hole2')
        profiles = []

        for pdb in self.find_files(self.pdb_dir, '_pos.pdb'):
            outfile = os.path.join('profiles', os.path.splitext(os.path.split(pdb)[1])[0] + '.dat')
            os.makedirs(os.path.dirname(outfile), exist_ok=True)
            log.info('Analyzing model {}'.format(pdb))
            hole2_script = 'coord {0} \nradius {1} \ncvect 0.0 0.0 1.0 \ncpoint 0.  0.  -15.\n'.format(pdb, hole2_rads)
            with open(outfile, 'w') as pore_file:
                pore_file.write(sp.check_output([hole2], input=hole2_script, universal_newlines=True))
                profiles.append(pore_file)
            log.info('Analysis done, results saved in {}'.format(outfile))

        return profiles

    def parse_profiles(self):
        """
        :return:
        """
        parsed_profiles = []
        min_rads = []

        for profile in self.find_files(self.pdb_dir, '_pos.dat'):
            with open(profile, 'r') as prof:
                parsed_profile = []
                for line in prof.readlines():
                    line = line.split()
                    if len(line) == 5 and line[4] in ['(sampled)', '(mid-point)']:
                        parsed_profile.append([float(col) for col in line[0:2]])
                    if len(line) == 5 and line[0] == 'Minimum':
                        min_rads.append(float(line[3]))
                parsed_profiles.append(np.array(parsed_profile))
        return [parsed_profiles, min_rads]

    def plot_profiles(self):
        """
        :return:
        """
        mmplt.plot_simple_multiple_numpy(self.parsed_profiles, "Z-axis [A]", "Radius [A]", self.labels,
                                         self.name + '_pore_profiles', sizex=3.75, sizey=3.0, ranges=True)

parser = ap.ArgumentParser()
parser.add_argument("-d", "--pdb_dirs", nargs='+', help="directories with pdb files to analyze")
parser.add_argument("-n", "--pdb_names", nargs='+', help="model names")
args = parser.parse_args()

lg.basicConfig(level=lg.INFO)
log = lg.getLogger('mm_analysis_pore')
models = PoreModel(args.pdb_dirs[0], args.pdb_names[0])
models.plot_profiles()
