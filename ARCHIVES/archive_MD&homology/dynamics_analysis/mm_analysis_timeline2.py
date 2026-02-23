# python 3 + pandas
# plots vmd rmsd tml files
# michaladammichalowski@gmail.com
#
#

import argparse
import mm_pandas_plot as mmpl
import mm_pandas_parse as mmpr

class Timeline:

    def __init__(self, in_files, out_file, property, segments, timestep=False, frequency=False):
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
        self.property = property
        self.segments = segments

        if self.timestep: self.timeframe = self.timestep*self.frequency*0.000001
        self.tml = self.read_tml()

    def read_tml(self):
        """
        :return: list of rmsd profiles
        """
        parser = mmpr.Parser

        par_tmls = parser.read_result_horizon(self.in_files, 9, [['a', 'b', 'c', 'd']])
        if self.timestep: parser.reindex_result(par_tmls, par_tmls.index*self.timeframe, 'time [ns]')

        print(par_tmls)

        return par_tmls

#    def plot_rmsd(self):
#        """
#        :return: plots parsed rdf profiles
#        """
#        ploter = mmpl.Ploter()
#        ploter.plot_single('multi', self.rmsds, self.out_file, [4.5, 2.5], y_label='RMSD [A]', ncols=3)

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--tml_files", nargs='+', help="rmsd files to plot")
parser.add_argument("-p", "--property", nargs='+', help="labels for data series")
parser.add_argument("-s", "--segments",  nargs='+', help="labels for data series")
parser.add_argument("-o", "--out_file", help="outfile name")
parser.add_argument("-t", "--timestep", type=float, help="timestep")
parser.add_argument("-f", "--frequency", type=int, help="dcd save frequency")

args = parser.parse_args()

tmls = Timeline(args.tml_files, args.out_file, args.property, args.segments, args.timestep, args.frequency)
