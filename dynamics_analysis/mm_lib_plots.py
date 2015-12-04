# python 3
# plot functions
# michaladammichalowski@gmail.com
# 25.11.15 - creation

import numpy as np
import matplotlib.pyplot as plt


def plot_simple(x, y, xlab, ylab, labe, out_file, color='#405952', fontsize=12, sizex=3.5, sizey=3.5):
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
    # ax.set_xlim([a, b])
    # ax.set_ylim([a, b])

    # LEGEND OPTIONS
    # handles, labels = ax.get_legend_handles_labels()
    # lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=fontsize)
    # lgd = ax.legend(handles, labels, loc='best', fontsize=fontsize)

    fig.set_size_inches(sizex, sizey)
    fig.savefig(out_file+".png", dpi=300, bbox_inches='tight')


def plot_simple_multiple(x, y, xlab, ylab, labe, out_file, fontsize=12, sizex=3.5, sizey=3.5):
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
    # ax.set_color_cycle(['#4D4D4D','#5DA5DA', '#60BD68', '#F17CB0', '#B2912F', '#B276B2', '#DECF3F', '#F15854'])
    ax.set_color_cycle(['#405952', '#9C9B7A', '#FFD393', '#FF974F', '#F54F29'])

    for serie in range(len(labe)):
        ax.plot(x[serie], y[serie], label=labe[serie], lw=3)

    ax.set_xlabel(xlab, fontsize=fontsize)
    ax.set_ylabel(ylab, fontsize=fontsize)
    ax.grid('on')
    ax.ticklabel_format(style='sci', scilimits=(-3, 4), axis='both')

    # RANGE OPTIONS:
    # ax.set_xlim([a, b])
    # ax.set_ylim([a, b])

    # LEGEND OPTIONS
    # handles, labels = ax.get_legend_handles_labels()
    # lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=fontsize)
    # lgd = ax.legend(handles, labels, loc='best', fontsize=fontsize)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=5, fontsize=fontsize)

    fig.set_size_inches(sizex, sizey)
    # fig.savefig(labe+".png", dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(out_file+".png", dpi=300, bbox_inches='tight')


def plot_histogram(data, xlab, ylab, labe, out_file, fontsize=12, sizex=3.5, sizey=3.5):
    """
    :param data:
    :param xlab:
    :param ylab:
    :param labe:
    :param out_file:
    :param fontsize:
    :param sizex:
    :param sizey:
    :return:
    """
    fig, ax = plt.subplots()
    # ax.set_color_cycle(['#4D4D4D','#5DA5DA', '#60BD68', '#F17CB0', '#B2912F', '#B276B2', '#DECF3F', '#F15854'])
    ax.set_color_cycle(['#405952', '#9C9B7A', '#FFD393', '#FF974F', '#F54F29'])

    for serie in range(len(labe)):
        n, bins, patches = plt.hist(data[serie], bins=150, label=labe[serie], normed=True, histtype='step', lw=3)

    ax.set_xlabel(xlab, fontsize=fontsize)
    ax.set_ylabel(ylab, fontsize=fontsize)
    ax.grid('on')
    ax.ticklabel_format(style='sci', scilimits=(-3, 4), axis='both')

    # RANGE OPTIONS:
    # ax.set_xlim([a, b])
    # ax.set_ylim([a, b])

    # LEGEND OPTIONS
    # handles, labels = ax.get_legend_handles_labels()
    # lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=fontsize)
    # lgd = ax.legend(handles, labels, loc='best', fontsize=fontsize)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=5, fontsize=fontsize)

    fig.set_size_inches(sizex, sizey)
    # fig.savefig(labe+".png", dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(out_file+".png", dpi=300, bbox_inches='tight')