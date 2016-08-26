import pandas as pd
import numpy as np


class Parser:

    def __init__(self):
        pass

    @staticmethod
    def read_result(fname, skip, header):
        """
        read space separated data file
        :param fname: file name
        :param skip: number of top rows to skip
        :param header: header (list of strings) for dataframe
        :return: dataframe with parsed data
        """
        return pd.read_csv(fname, delim_whitespace=True, skiprows=skip, names=header, index_col=0)
