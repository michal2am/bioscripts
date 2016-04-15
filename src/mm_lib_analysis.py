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


def filter_out(values, index, how, treshold):
    """
    :param values: line-like values iterable
    :param index: property to check
    :param how: up (remove higher), bottom (remove lower), out (remove outside range)
    :param treshold: limit(s) to remove
    :return:
    """
    if how == 'up':
        compare = lambda x, trsh: True if x < trsh else False
    elif how == 'bottom':
        compare = lambda x, trsh: True if x > trsh else False
    elif how == 'out':
        compare = lambda x, trsh: True if (x > trsh[0] and x < trsh[1]) else False

    return [element[:] for element in values if compare(float(element[index]), treshold)]

