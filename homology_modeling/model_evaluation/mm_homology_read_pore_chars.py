# python 3
# parser of caver 3.0 results
# michaladammichalowski@gmail.com
# 26.12.15 - creation

import csv
import argparse


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
        if tunnel_lenght > 80:
            print("Tunnel number:{0} length: {1}".format(tunnel_chars["Tunnel"][tunnel_chars["Length"].index(tunnel_lenght)], tunnel_lenght))
            long_tunnels.append(tunnel_chars["Tunnel"][tunnel_chars["Length"].index(tunnel_lenght)])
        if tunnel_shape > 2.0:
            print("Tunnel number:{0} curvature: {1}".format(tunnel_chars["Tunnel"][tunnel_chars["Curvature"].index(tunnel_shape)], tunnel_shape))
            stra_tunnels.append(tunnel_chars["Tunnel"][tunnel_chars["Curvature"].index(tunnel_shape)])

    for tunnel in long_tunnels:
        if tunnel in stra_tunnels:
            print("Good tunnel: {}".format(tunnel))

parser = argparse.ArgumentParser()
parser.add_argument("-tc", "--tunnel_chars", help="caver 3.0 tunnel characteristics in .csv format")
args = parser.parse_args()

tunnel_chars = (parse_tunnel_characteristics(args.tunnel_chars))
find_main_tunnels(tunnel_chars)





