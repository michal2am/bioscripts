import pandas as pd
import numpy as np


class Parser:

    def __init__(self):
        pass

    @staticmethod
    def read_result_vert(fname, skip, header):
        """
        read space separated data file
        :param fname: file name or names
        :param skip: number of top rows to skip
        :param header: header (list of strings) for dataframe
        :return: dataframe with parsed data
        """
        if isinstance(fname, str):
            print('in string')
            return pd.read_csv(fname, delim_whitespace=True, skiprows=skip, header=None, names=header, index_col=0)

        else:
            print(fname, header)
            return pd.concat([Parser.read_result_vert(fn, skip, [hd]) for (hd, fn) in zip(header, fname)], axis=1)

    @staticmethod
    def read_result_horizon(fname, skip, header):
        """

        :param fname:
        :param skip:
        :param header:
        :return:
        """
        if isinstance(fname, str):
            a = pd.read_csv(fname, delim_whitespace=True, skiprows=skip, header=None, names=header, index_col=0)
            print('in string')
            return pd.read_csv(fname, delim_whitespace=True, skiprows=skip, header=None, names=header, index_col=0)

        else:
            print(fname, header)
            return pd.concat([Parser.read_result_horizon(fn, skip, [hd]) for (hd, fn) in zip(header, fname)], axis=1)

    @staticmethod
    def reindex_result(dframe, index, header):
        """

        :param dframe:
        :param func:
        :param header:
        :return:
        """
        dframe.set_index(index, inplace=True, drop=True)
        dframe.index.name = header

        return dframe
