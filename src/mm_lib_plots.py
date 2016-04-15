# python 3
# plot functions
# michaladammichalowski@gmail.com
# 25.11.15 - creation

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.collections as cls
import matplotlib.ticker as tck


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

    # SCALE OPTIONS:
    # ax.set_xscale('log')
    # ax.set_yscale('log')

    # RANGE OPTIONS:
    # ax.set_xlim([a, b])
    # ax.set_ylim([a, b])

    # LEGEND OPTIONS
    # handles, labels = ax.get_legend_handles_labels()
    # lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=fontsize)
    # lgd = ax.legend(handles, labels, loc='best', fontsize=fontsize)

    fig.set_size_inches(sizex, sizey)
    fig.savefig(out_file+".png", dpi=300, bbox_inches='tight')


def plot_simple_multiple(x, y, xlab, ylab, labe, out_file, ylimit=False, fontsize=12, sizex=3.5, sizey=3.5, linestyle='-', marker='None'):
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
        ax.plot(x[serie], y[serie], label=labe[serie], lw=2, linestyle=linestyle, marker=marker)

    ax.set_xlabel(xlab, fontsize=fontsize)
    ax.set_ylabel(ylab, fontsize=fontsize)
    ax.grid('on')
    ax.ticklabel_format(style='sci', scilimits=(-3, 4), axis='both')

    # RANGE OPTIONS:
    if ylimit:
        ax.set_ylim(ylimit[0], ylimit[1])
    # ax.set_xlim([a, b])

    # LEGEND OPTIONS
    # handles, labels = ax.get_legend_handles_labels()
    # lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=fontsize)
    # lgd = ax.legend(handles, labels, loc='best', fontsize=fontsize)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=fontsize)

    fig.set_size_inches(sizex, sizey)
    # fig.savefig(labe+".png", dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(out_file+".png", dpi=300, bbox_inches='tight')

def plot_simple_multiple_numpy(data, xlab, ylab, labe, out_file, ranges=False, xlimit=False, ylimit=False, back=False, fontsize=6, sizex=3.5, sizey=3.5, linestyle='-', marker='None'):
    """
    :param data: list of numpy arrays with rows for single data points
    :param xlab: x-axis label
    :param ylab: y-axis label
    :param labe: data label
    :param color: data color
    :param fontsize:
    :param ranges: turn on for colormap scaling instead of "random" colors
    :return:
    """
    fig, ax = plt.subplots()
    colors = ['#9C9B7A', '#405952', '#703030', '#2F343B','#FFD393', '#FF974F', '#F54F29', '#B38F73', '#3A4012', '#C4D8F2', '#5F6F82',
              '#3E4048', '#C1DBBD', '#6B85B5', '#30507E']
    if ranges:
        color_map = mpl.cm.get_cmap('inferno')
        norm = mpl.colors.Normalize(vmin=0, vmax=len(labe))

    if back:
        c1 = cls.BrokenBarHCollection(back[0], back[1], facecolor='orange', alpha=0.3)
        ax.add_collection(c1)

    for serie in range(len(labe)):
        color = color_map(norm(serie)) if ranges else colors[serie]
        ax.plot(data[serie][:, 0], data[serie][:, 1], label=labe[serie], lw=0.5, linestyle=linestyle,
                marker=marker, color=color)
        for adds in range(data[serie].shape[1] - 2):
            ax.plot(data[serie][:, 0], data[serie][:, adds + 2], label='_nolegend_', lw=0.5, linestyle=linestyle,
                    marker=marker, color=colors[serie])

    ax.set_xlabel(xlab, fontsize=fontsize)
    ax.set_ylabel(ylab, fontsize=fontsize)
    ax.grid('on')

    ### TICK CONFIGURATION ###

    ax.tick_params(labelsize=fontsize)
    ax.ticklabel_format(style='sci', scilimits=(-3, 4), axis='both')

    # defaul tick values
    ax.yaxis.set_major_locator(tck.MultipleLocator(2))
    ax.yaxis.set_minor_locator(tck.MultipleLocator(0.5))
    ax.xaxis.set_major_locator(tck.MultipleLocator(50))
    ax.xaxis.set_minor_locator(tck.MultipleLocator(12.5))

    if ylimit:
        ax.set_ylim(ylimit[0], ylimit[1])
    if xlimit:
        ax.set_xlim(xlimit[0], xlimit[2])

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=2, fontsize=fontsize)

    fig.set_size_inches(sizex, sizey)
    fig.savefig(out_file+".png", dpi=300, bbox_inches='tight')


def plot_histogram(data, xlab, ylab, labe, out_file, fontsize=12, sizex=3.5, sizey=3.5, ranges=False):
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

    colors = ['#9C9B7A', '#405952', '#703030', '#2F343B','#FFD393', '#FF974F', '#F54F29', '#B38F73', '#3A4012', '#C4D8F2', '#5F6F82',
              '#3E4048', '#C1DBBD', '#6B85B5', '#30507E']

    if ranges:
        color_map = mpl.cm.get_cmap('seismic')
        norm = mpl.colors.Normalize(vmin=0, vmax=len(labe))

    for serie in range(len(labe)):
        color = color_map(norm(serie)) if ranges else colors[serie]
        n, bins, patches = plt.hist(data[serie], orientation='horizontal', bins=20, label=labe[serie], normed=True, histtype='step', lw=0.7, color=color)

    ### LABEL CONFIGURATION ###

    ax.set_xlabel(xlab, fontsize=fontsize)
    ax.set_ylabel(ylab, fontsize=fontsize)

    ax.grid('on')

    ### TICK CONFIGURATION ###

    ax.tick_params(labelsize=fontsize)
    ax.ticklabel_format(style='sci', scilimits=(-3, 4), axis='both')

    # defaul tick values
    ax.yaxis.set_major_locator(tck.MultipleLocator(10))
    ax.yaxis.set_minor_locator(tck.MultipleLocator(2.5))
    ax.xaxis.set_major_locator(tck.MultipleLocator(0.01))
    ax.xaxis.set_minor_locator(tck.MultipleLocator(0.0025))

    # custom x tick values
    # xticks = [int(tick) for tick in np.arange(min(data[0]), max(data[0]), 20)]
    # ax.set_xticks(xticks)
    # ax.set_xticklabels(xticks)

    ### RANGE OPTIONS ###

    # ax.set_xlim([a, b])
    # ax.set_ylim([a, b])

    ### LEGEND OPTIONS ###

    # handles, labels = ax.get_legend_handles_labels()
    # lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=fontsize)
    # lgd = ax.legend(handles, labels, loc='best', fontsize=fontsize)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.4), ncol=2, fontsize=fontsize)

    fig.set_size_inches(sizex, sizey)
    # fig.savefig(labe+".png", dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(out_file+".png", dpi=300, bbox_inches='tight')


def plot_ticker(data, xlab, ylab, out_file, color='#405952', fontsize=10, sizex=3.5, sizey=3.5):
    """
    :param data:
    :param xlab:
    :param ylab:
    :param out_file:
    :param color:
    :param fontsize:
    :param sizex:
    :param sizey:
    :return:
    """

    fig, ax = plt.subplots()
    x = list(range(len(data)))
    y = []
    labels = []
    for key in data:
        labels.append(key)
        y.append(data[key])

    ax.bar(x, y, align='center', width=1, color=color)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation='vertical', fontsize=fontsize)
    ax.set_xlabel(xlab, fontsize=fontsize)
    ax.set_ylabel(ylab, fontsize=fontsize)
    ax.grid('on')

    # SCALE OPTIONS:
    # ax.set_xscale('log')
    # ax.set_yscale('log')

    # RANGE OPTIONS:
    ax.set_xlim([-1, len(x)])
    # ax.set_ylim([a, b])

    # LEGEND OPTIONS
    # handles, labels = ax.get_legend_handles_labels()
    # lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=fontsize)
    # lgd = ax.legend(handles, labels, loc='best', fontsize=fontsize)

    fig.set_size_inches(sizex, sizey)
    fig.savefig(out_file+".png", dpi=300, bbox_inches='tight')
