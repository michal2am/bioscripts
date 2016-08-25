import pandas as pd
import numpy as np


class TimePositions:

    def __init__(self, filname):
        self.filename = filname
        self.positions = self.read_result(self.filename, 0, ['xr', 'yr', 'zr', 'x', 'y', 'z'])


    @staticmethod
    def euclidean(values):
        print(values)
        print(type(values))
        print(values.iloc[0:3])
        print(values.iloc[3:6])
        print(values.iloc[0:3].subtract(values.iloc[3:6]))

        return np.linalg.norm(values.iloc[0:2] - values.iloc[3:5])

    @staticmethod
    def read_result(fname, skip, header):
        """
        read space separated data file
        :param fname: file name
        :param skip: number of top rows to skip
        :param header: header (list of strings) for dataframe
        :return: dataframe with parsed data
        """

        data = pd.read_csv(fname, delim_whitespace=True, skiprows=skip, names=header, index_col=0)
        return data

filename = 'z_chlorides.dat'
results = TimePositions(filename)
#print(results.positions)
print(results.positions.apply(results.euclidean, axis=1))
