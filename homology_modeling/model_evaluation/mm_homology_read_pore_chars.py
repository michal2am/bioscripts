# python 3
# parser of caver 3.0 results
# michaladammichalowski@gmail.com
# 26.12.15 - creation

import csv
import argparse
import numpy as np
import mm_lib_plots as mmplt


def parse_tunnel_characteristics(csv_tunnel_characteristics):
    """
    :param csv_tunnel_characteristics:
    :return:
    """
    tunnel_chars_keys = ["Snapshot", "Tunnel cluster", "Tunnel", "Throughput", "Cost", "Bottleneck radius", \
                         "Bottleneck R error bound", "Length", "Curvature"]
    tunnel_chars = {key: [] for key in tunnel_chars_keys}
    dict_reader = csv.DictReader(open(csv_tunnel_characteristics, 'r'), fieldnames=tunnel_chars_keys, \
                                 delimiter=',', quotechar='"')
    next(dict_reader)
    for row in dict_reader:
        for key in row:
            if key == "Snapshot":
                tunnel_chars[key].append(row[key])
            elif key in ["Tunnel cluster", "Tunnel"]:
                tunnel_chars[key].append(int(row[key]))
            else:
                tunnel_chars[key].append(float(row[key]))

    return tunnel_chars


def find_main_tunnels(tunnel_chars):
    """
    :param tunnel_chars:
    :return:
    """
    long_tunnels = []
    stra_tunnels = []

    for tunnel_lenght, tunnel_shape in zip(tunnel_chars["Length"], tunnel_chars["Curvature"]):
        if tunnel_lenght > 80.0:
            print("Tunnel number:{0} length: {1}".format(tunnel_chars["Tunnel"][tunnel_chars["Length"].index(tunnel_lenght)], tunnel_lenght))
            long_tunnels.append(tunnel_chars["Tunnel"][tunnel_chars["Length"].index(tunnel_lenght)])
        if tunnel_shape < 1.2:
            print("Tunnel number:{0} curvature: {1}".format(tunnel_chars["Tunnel"][tunnel_chars["Curvature"].index(tunnel_shape)], tunnel_shape))
            stra_tunnels.append(tunnel_chars["Tunnel"][tunnel_chars["Curvature"].index(tunnel_shape)])

    good_tunels = []

    for tunnel in long_tunnels:
        if tunnel in stra_tunnels:
            good_tunels.append(tunnel)
            print("Good tunnel: {}".format(tunnel))

    return good_tunels


def parse_selected_profiles(csv_tunnel_profiles, selected_tunnels):
    """
    :param selected_tunnels:
    :return:
    """
    profiles = [line.split(sep=',') for line in open(csv_tunnel_profiles, 'r').readlines()[1:-1]]
    selected_dist = []
    selected_rads = []

    for tunnel in selected_tunnels:
        for line in profiles:
            if int(line[1]) == tunnel and line[12].strip() == 'R':
                selected_rads.append([float(rad) for rad in line[13:]])
            if int(line[1]) == tunnel and line[12].strip() == 'Z':
                selected_dist.append([float(rad) for rad in line[13:]])

    return [selected_dist, selected_rads]


def plot_profile(parseds_profs):
    """
    :param parseds_profs: list of parsed rdf profiles
    :return: plots parsed rdf profiles (no return)
    """
    mmplt.plot_simple_multiple(parseds_profs[1], parseds_profs[0], "distance [A]", "density [rel]", ['19', '20'], 'tmp.png')


parser = argparse.ArgumentParser()
parser.add_argument("-tc", "--tunnel_chars", help="caver 3.0 tunnel characteristics in .csv format")
parser.add_argument("-tp", "--tunnel_profiles", help="caver 3.0 tunnel profiles in .csv format")
args = parser.parse_args()

# tunnel_chars = (parse_tunnel_characteristics(args.tunnel_chars))
# find_main_tunnels(tunnel_chars)
parseds_profs = parse_selected_profiles(args.tunnel_profiles, [19, 20])
plot_profile(parseds_profs)




