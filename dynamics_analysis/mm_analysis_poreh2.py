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
import mm_lib_analysis as mmanl
import pandas as pd
import mm_pandas_parse as mmpr
import mm_pandas_plot as mmpl

# hole2 fixed file paths
hole2 = '/opt/hole2/exe/hole'
hole2_rads = '/opt/hole2/rad/simple.rad'
positioner = '/home/mm/Pycharm/bioscripts/homology_modeling/model_positioning/complete_positioner.sh'

class PoreModel:

    def __init__(self, pdb_dir, name):

        # modules stuff
        self.parser = mmpr.Parser()
        self.ploter = mmpl.Ploter()

        self.pdb_dir = pdb_dir
        self.name = name

        # self.set_position()
        # self.get_profile()
        self.labels = self.set_names()
        self.parsed_profiles, self.min_Rads = self.parse_profiles()


    @staticmethod
    def find_files_rec(directory, extension):
        """
        recursive
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

    @staticmethod
    def find_files_cur(directory, extension):
        """
        one level directory
        :param directory: directory to search in
        :param extension: files to look for
        :return: list of relative paths to files sorted by modification date
        """
        file_list = os.listdir(directory)
        file_list.sort(key=lambda x: os.stat(os.path.join(directory, x)).st_mtime)
        file_list = [file for file in file_list if file.endswith(extension)]
        return file_list

    def set_names(self):
        """
        :return: list of pdb file names without extension
        """
        positioned = [os.path.splitext(os.path.split(pdb)[1])[0] for pdb in self.find_files_cur(self.pdb_dir + '/profiles', '.dat')]
        return [name.split('_')[0] + ' ' + name.split('_')[1] for name in positioned]
        # return positioned

    def set_position(self):
        """
        calls complete_positioner to align and center models, pdb save in "positioned" dir
        """
        log.info('Searching in dirs: {} for models {}'.format(self.pdb_dir, self.name))
        log.info('Positioning models according to principal axes and geometry center')
        sp.check_call([positioner, self.pdb_dir])

    def get_profile(self):
        """
        :return: list of pore profiles files in "profiles" dir
        """
        log.info('Reading positioned models from {}'.format(self.pdb_dir))
        log.info('Esitmating pore profiles with Hole2')
        profiles = []

        for pdb in self.find_files_cur(self.pdb_dir + '/positioned', '_pos.pdb'):
            outfile = os.path.join('profiles', os.path.splitext(os.path.split(pdb)[1])[0] + '.dat')
            os.makedirs(os.path.dirname(outfile), exist_ok=True)
            log.info('Analyzing model {}'.format(pdb))
            pdb = os.path.join('positioned', pdb)
            hole2_script = 'coord {0} \nradius {1} \ncvect 0.0 0.0 1.0 \ncpoint 0.  0.  -15.\n'.format(pdb, hole2_rads)
            with open(outfile, 'w') as pore_file:
                pore_file.write(sp.check_output([hole2], input=hole2_script, universal_newlines=True))
                profiles.append(pore_file)
            log.info('Analysis done, results saved in {}'.format(outfile))

        return profiles

    def parse_profiles(self):
        """
        :return: list of parsed profiles and list of minimum radii
        """
        parsed_profiles = []
        min_rads = []

        for profile, name in zip(self.find_files_cur(self.pdb_dir + '/profiles', '_pos.dat'), self.labels):
            with open(os.path.join('profiles', profile), 'r') as prof:
                parsed_profile = []
                for line in prof.readlines():
                    line = line.split()
                    if len(line) == 5 and line[4] in ['(sampled)', '(mid-point)']:
                        parsed_profile.append([float(col) for col in line[1::-1]])
                    if len(line) == 5 and line[0] == 'Minimum':
                        min_rads.append(float(line[3]))
                parsed_profile = np.array(mmanl.filter_out(parsed_profile, 1, 'out', [-55, -10]))
                parsed_profile = pd.DataFrame(index=parsed_profile[:, 0], data=parsed_profile[:, 1], columns=[name])
                parsed_profile.index.name = 'radius [A]'
                parsed_profiles.append(parsed_profile)
        print(parsed_profiles)
        return [parsed_profiles, min_rads]

    def plot_profiles(self):
        """
        :return:
        """
        #mmplt.plot_simple_multiple_numpy(self.parsed_profiles, "Radius [A]", "Z-axis [A]", self.labels,
        #                                 self.name + '_pore_profiles', sizex=3.75, sizey=3.0, ranges=True,
        #                                 xlimit=[-0.5, 6.5])
        self.ploter.plot_multiple(self.parsed_profiles, 'pore profiles')

parser = ap.ArgumentParser()
parser.add_argument("-d", "--pdb_dirs", nargs='+', help="directories with pdb files to analyze")
parser.add_argument("-n", "--pdb_names", nargs='+', help="model names")
args = parser.parse_args()

lg.basicConfig(level=lg.INFO)
log = lg.getLogger('mm_analysis_pore')
models = PoreModel(args.pdb_dirs[0], args.pdb_names[0])
models.plot_profiles()
