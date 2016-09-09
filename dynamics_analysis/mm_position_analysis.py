# python ~/Pycharm/bioscripts/dynamics_analysis/mm_position_analysis.py -f xy_upring.dat -s 0.01 -l 'time [ns]' 'x1' 'x2' 'x3' 'x4' 'x5' 'y1' 'y2' 'y3' 'y4' 'y5' -t shoelace
# python ~/Pycharm/bioscripts/dynamics_analysis/mm_position_analysis.py -f z_chlorides.dat -s 0.01 -l 'time [ns]' 'xr' 'yr' 'zr' 'x' 'y' 'z' -t euclidean


import mm_pandas_parse as mmpr
import mm_pandas_plot as mmpl
import numpy as np
import argparse as arp
import pandas as pd


class TimePositions:

    def __init__(self, filname, step, columns, task, options, **kwargs):
        """
        for vmd .dat files of step + coords format parsing and plotting
        :param filname: string
        :param step: float, ns for step
        :param columns: list, columns names
        :param task: string, selected action
        :param options: list, optional
        :param annotations: TODO
        """

        # modules stuff
        self.parser = mmpr.Parser()
        self.ploter = mmpl.Ploter()

        # kwargs handling
        self.annot_points, self.annot_texts = kwargs['annotations'][0], kwargs['annotations'][1]
        self.bounds = [float(b) for b in kwargs['bounds']] if kwargs['bounds'] is not None else None

        # basic parsing
        self.filename, self.step, self.options = filname, step, options
        self.positions = self.parser.read_result(self.filename, 0, columns)
        if self.bounds is not None:
            self.positions = self.positions.loc[self.bounds[0]:self.bounds[1]]
            label = self.positions.index.name
            self.positions.reset_index(drop=True, inplace=True)
            self.positions.index.name = label
        self.tons()

        if task == 'euclidean':
            self.euclidean()
            self.plot_euclidean()
            self.plot_n()
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
        self.ploter.plot_single('single', self.positions.loc[:, 'area [A^2]'].to_frame(), 'upring area [A^2]', (4, 3),
                                annot_points=self.annot_points, annot_texts=self.annot_texts, bounds=self.bounds)

    def plot_euclidean(self):
        """
        plots euclidean distance time plot
        """
        self.ploter.plot_single('single', self.positions.loc[:, 'distance [A]'].to_frame(), 'chloride positions', (4, 3),
                                annot_points=self.annot_points, annot_texts=self.annot_texts, bounds=self.bounds)

    def plot_both(self):
        self.ploter.plot_multiple('a', [self.positions.loc[:, 'area [A^2]'].to_frame(), self.positions.loc[:, 'distance [A]'].to_frame()], (4, 6))

    def plot_n(self):
        self.ploter.plot_single('single', self.positions.loc[:, 'n'].to_frame(), 'chloride ions within 4A', (4, 3),
                                annot_points=self.annot_points, annot_texts=self.annot_texts, bounds=self.bounds)


parser = arp.ArgumentParser()
parser.add_argument("-f", "--file", help="file name")
parser.add_argument("-s", "--step", type=float, help="file time step")
parser.add_argument("-l", "--labels", nargs="+", help="file data labels")
parser.add_argument("-t", "--task", help="action to perform")
parser.add_argument("-ap", "--annot_point", nargs="+", help="plot annotations: points")
parser.add_argument("-at", "--annot_text", nargs="+", help="plot annotations: titles")
parser.add_argument("-b", "--bounds", nargs="+", help="start and stop timeframe to plot")
parser.add_argument("-o", "--options", nargs="+", help="additional options")
args = parser.parse_args()

results = TimePositions(args.file, args.step, args.labels, args.task, args.options,
                        annotations=(args.annot_point, args.annot_text), bounds=args.bounds)
