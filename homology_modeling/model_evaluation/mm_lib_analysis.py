# python 3
# simple analysis functions
# michaladammichalowski@gmail.com
# 24.02.16 - creation

import numpy as np
import matplotlib.pyplot as plt


def get_sort(values, index):
    """
    :param values: line-like values iterable
    :param index: property to compare
    :return: sorted line-like values iterable with position at [0]
    """
    return [[i+1] + element[:] for i, element in enumerate(sorted(values, key=lambda x: x[index]))]

