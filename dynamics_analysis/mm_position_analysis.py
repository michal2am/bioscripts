import mm_pandas_parse as mmpr
import mm_pandas_plot as mmpl
import numpy as np
import pandas as pd


class TimePositions:

    def __init__(self, filname, step, columns):
        self.parser = mmpr.Parser()
        self.ploter = mmpl.Ploter()

        self.filename, self.step = filname, step
        self.positions = self.parser.read_result(self.filename, 0, columns)
        self.tons()
        self.euclidean()
        #self.shoelace()

    def tons(self):
        self.positions.index = self.positions.index * self.step

    def euclidean(self):
        """
        adds dist column with euclidean distance
        for dataframe ['time [ns]', 'xr', 'yr', 'zr', 'x', 'y', 'z']
        :return: none
        """
        self.positions['distance [A]'] = np.sqrt((self.positions.x - self.positions.xr)**2 +\
                                         (self.positions.y - self.positions.yr)**2 +\
                                         (self.positions.z - self.positions.zr)**2)

    def shoelace(self):
        def poly_area(x, y):
            return 0.5*np.abs(np.dot(x, np.roll(y, 1))-np.dot(y, np.roll(x, 1)))
        self.positions['area [A^2]'] = [poly_area(timeframe[0:5], timeframe[5:10]) for timeframe in self.positions.values]


    def plot_shoelace(self):
        self.ploter.plot_vertical_result([self.positions[['area [A^2]']]], ['line'], 'upring area [A^2]')

    def plot_euclidean(self):
        """
        plots euclidean distance time plot and histogram
        :return:
        """
        self.ploter.plot_single(self.positions[['distance [A]']], 'line', 'chloride positions')


filename = 'z_chlorides.dat'
step = 0.01
columns = ['time [ns]', 'xr', 'yr', 'zr', 'x', 'y', 'z']
results = TimePositions(filename, step, columns)
results.plot_euclidean()
'''
filename = 'xy_upring.dat'
step = 0.01
columns = ['time [ns]', 'x1', 'x2', 'x3', 'x4', 'x5', 'y1', 'y2', 'y3', 'y4', 'y5']
results = TimePositions(filename, step, columns)
print(results.positions)
results.plot_shoelace()
'''