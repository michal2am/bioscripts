# python 3
# plot functions
# michaladammichalowski@gmail.com
# 25.11.15 - creation

import numpy as np
import matplotlib.pyplot as plt


def plot_simple(x, y, xlab, ylab, labe, color='#f69855', fontsize=18, sizex=3.5, sizey=3.5):
    """
    :param x: x-axis data
    :param y: y-axis data
    :param xlab: x-axis label
    :param ylab: y-axis label
    :param labe: data label
    :param color: data color
    :param fontsize:
    :return:
    """
    fig, ax = plt.subplots()
    ax.plot(x, y, label=labe, lw=3, color=color)
    ax.set_xlabel(xlab, fontsize=fontsize)
    ax.set_ylabel(ylab, fontsize=fontsize)
    ax.grid('on')
    ax.ticklabel_format(style='sci', scilimits=(-3, 4), axis='both')

    # RANGE OPTIONS:
    # ax.set_xlim([0, 120])
    # ax.set_ylim([0.75, 1.02])

    # LEGEND OPTIONS
    # handles, labels = ax.get_legend_handles_labels()
    # lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=10)
    # lgd = ax.legend(handles, labels, loc='best')

    fig.set_size_inches(sizex, sizey)
    # fig.savefig(labe+".png", dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(labe+".png", dpi=300, bbox_inches='tight')