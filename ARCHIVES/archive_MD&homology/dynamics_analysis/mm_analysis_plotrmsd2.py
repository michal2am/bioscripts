# python 3 + pandas
# plots vmd rmsd dat files
# michaladammichalowski@gmail.com
# 21.12.15 - creation
# 10.03.15 - refactor
# 15.11.16 - refactor
# python3 ~/pycharm_projects/bioscripts/dynamics_analysis/mm_analysis_plotrmsd2.py --rmsd_files ALL.dat  A_a.dat B_a.dat C_a.dat D_a.dat E_a.dat --labels all A B C D E --out_file rmsd_time_plot --timestep 2.0 --frequency 25000  --ranges 0 5

import argparse
import mm_pandas_plot as mmpl
import mm_pandas_parse as mmpr

class RMSD:

    def __init__(self, in_files, out_file, labels, timestep=False, frequency=False):
        """
        :param in_files: rmsd files to read
        :param out_file: file to save plot
        :param timestep: optional timestep for real time axis
        :param frequency: optional dcd frequency for real time axis
        :param labels: data series labels
        :return:
        """
        self.in_files = in_files
        self.out_file = out_file
        self.timestep = timestep
        self.frequency = frequency
        self.labels = labels

        if self.timestep: self.timeframe = self.timestep*self.frequency*0.000001
        self.rmsds = self.read_rmsd()

    def read_rmsd(self):
        """
        :return: list of rmsd profiles
        """
        parser = mmpr.Parser

        par_rmsds = parser.read_result_vert(self.in_files, 2, self.labels)
        if self.timestep: parser.reindex_result(par_rmsds, par_rmsds.index*self.timeframe, 'time [ns]')

        return par_rmsds

    def plot_rmsd(self):
        """
        :return: plots parsed rdf profiles
        """
        ploter = mmpl.Ploter()
        ploter.plot_single('multi', self.rmsds, self.out_file, [5.0, 2.5], y_label='RMSD [A]', ncols=3)

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--rmsd_files", nargs='+', help="rmsd files to plot")
parser.add_argument("-l", "--labels", nargs='+', help="labels for data series")
parser.add_argument("-o", "--out_file", help="outfile name")
parser.add_argument("-t", "--timestep", type=float, help="timestep")
parser.add_argument("-f", "--frequency", type=int, help="dcd save frequency")
parser.add_argument("-a", "--ranges", nargs='+', type=float, help="RMSD plot range")

args = parser.parse_args()

rmsds = RMSD(args.rmsd_files, args.out_file, args.labels, args.timestep, args.frequency)
rmsds.plot_rmsd()
