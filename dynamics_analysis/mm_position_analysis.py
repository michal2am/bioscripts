# python ~/Pycharm/bioscripts/dynamics_analysis/mm_position_analysis.py -f xy_upring.dat -s 0.01 -l 'time [ns]' 'x1' 'x2' 'x3' 'x4' 'x5' 'y1' 'y2' 'y3' 'y4' 'y5' -t shoelace
# python ~/Pycharm/bioscripts/dynamics_analysis/mm_position_analysis.py -f z_chlorides.dat -s 0.01 -l 'time [ns]' 'xr' 'yr' 'zr' 'x' 'y' 'z' -t euclidean


import mm_pandas_parse as mmpr
import mm_pandas_plot as mmpl
import numpy as np
import argparse as arp


class TimePositions:

    def __init__(self, filname, step, columns, task, options):
        """
        for vmd .dat files of step + coords format parsing and plotting
        :param filname: string
        :param step: float, ns for step
        :param columns: list, columns names
        :param task: string, selected action
        :param options: list, optional
        """
        self.parser = mmpr.Parser()
        self.ploter = mmpl.Ploter()

        self.filename, self.step, self.options = filname, step, options
        self.positions = self.parser.read_result(self.filename, 0, columns)
        self.tons()

        if task == 'euclidean':
            self.euclidean()
            print(self.positions.loc[2.17].name)
            self.plot_euclidean()
        if task == 'shoelace':
            self.shoelace()
            self.plot_shoelace()

    def tons(self):
        """
        converts index to ns timescale
        """
        self.positions.index = self.positions.index * self.step

    def euclidean(self):
        """
        adds dist column with euclidean distance
        for dataframe ['time [ns]', 'xr', 'yr', 'zr', 'x', 'y', 'z']
        """
        self.positions['distance [A]'] = np.sqrt((self.positions.x - self.positions.xr)**2 +\
                                         (self.positions.y - self.positions.yr)**2 +\
                                         (self.positions.z - self.positions.zr)**2)

    def shoelace(self):
        """
        adds area column with shoelace area
        for dataframe ['time [ns]', 'x1', 'x2', ... 'y1', 'y2']
        """
        def poly_area(x, y):
            return 0.5*np.abs(np.dot(x, np.roll(y, 1))-np.dot(y, np.roll(x, 1)))
        point_num = int(self.positions.shape[1]/2)
        self.positions['area [A^2]'] = [poly_area(timeframe[0:point_num], timeframe[point_num:]) for timeframe in self.positions.values]

    def plot_shoelace(self):
        """
        plots shoelace area time plot
        """
        self.ploter.plot_single(self.positions.loc[:, 'area [A^2]'].to_frame(), 'upring area [A^2]')

    def plot_euclidean(self):
        """
        plots euclidean distance time plot
        """
        self.ploter.plot_single(self.positions.loc[:, 'distance [A]'].to_frame(), 'chloride positions')


parser = arp.ArgumentParser()
parser.add_argument("-f", "--file", help="file name")
parser.add_argument("-s", "--step", type=float, help="file time step")
parser.add_argument("-l", "--labels", nargs="+", help="file data labels")
parser.add_argument("-t", "--task", help="action to perform")
parser.add_argument("-o", "--options", nargs="+", help="additional options")
args = parser.parse_args()

results = TimePositions(args.file, args.step, args.labels, args.task, args.options)
