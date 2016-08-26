import mm_pandas_parse as mmpr
import mm_pandas_plot as mmpl
import numpy as np
import pandas as pd


class TimePositions:

    def __init__(self, filname, step):
        self.parser = mmpr.Parser()
        self.ploter = mmpl.Ploter()

        self.filename, self.step = filname, step
        self.positions = self.parser.read_result(self.filename, 0, ['time [ns]', 'xr', 'yr', 'zr', 'x', 'y', 'z'])
        self.tons()
        self.euclidean()

    def tons(self):
        self.positions.index = self.positions.index * self.step

    def euclidean(self):
        """
        adds dist column with euclidean distance
        :return: none
        """
        self.positions['distance [A]'] = np.sqrt((self.positions.x - self.positions.xr)**2 +\
                                         (self.positions.y - self.positions.yr)**2 +\
                                         (self.positions.z - self.positions.zr)**2)

    def plot(self):
        """
        plots euclidean distance time plot and histogram
        :return:
        """
        self.ploter.plot_vertical_result([self.positions[['distance [A]']], self.positions[['distance [A]']]], ['line', 'hist'], 'chloride positions')


filename = 'z_chlorides.dat'
step = 0.01
results = TimePositions(filename, step)
results.plot()
