# python 3
# parser of caver 3.0 results
# michaladammichalowski@gmail.com
# 26.12.15 - creation
# EXAMPLE CALL1: python3 mm_homology_read_pore_chars.py --tunnel_chars mod56/analysis/tunnel_characteristics.csv --tunnel_profiles mod56/analysis/tunnel_profiles.csv --outfile profiles --range 30 80
# EXAMPLE CALL2: python3 mm_homology_read_pore_chars.py --tunnel_chars mod56/analysis/tunnel_characteristics.csv --tunnel_profiles mod56/analysis/tunnel_profiles.csv --outfile profiles --range 30 80 --pdbs 3JAE_pos.pdb --profiles 19 --labels "GlyR + gly"
# NOTES: assuming pore profile starting from bottom (negative Z)

import csv
import argparse
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
        if tunnel_lenght > 70.0:
            print("Tunnel number:{0} length: {1}".format(tunnel_chars["Tunnel"][tunnel_chars["Length"].index(tunnel_lenght)], tunnel_lenght))
            long_tunnels.append(tunnel_chars["Tunnel"][tunnel_chars["Length"].index(tunnel_lenght)])
        if tunnel_shape < 1.2:
            print("Tunnel number:{0} curvature: {1}".format(tunnel_chars["Tunnel"][tunnel_chars["Curvature"].index(tunnel_shape)], tunnel_shape))
            stra_tunnels.append(tunnel_chars["Tunnel"][tunnel_chars["Curvature"].index(tunnel_shape)])

    good_tunels = []
    snapshots = []

    for tunnel in long_tunnels:
        if tunnel in stra_tunnels:
            snapshot = (tunnel_chars["Snapshot"][tunnel_chars["Tunnel"].index(tunnel)])
            snapshots.append(snapshot)
            good_tunels.append(tunnel)
            print("{0} good tunnel: {1}".format(snapshot, tunnel))

    return [snapshots, good_tunels]


def parse_selected_profiles(csv_tunnel_profiles, snapshots, selected_tunnels, range):
    """
    :param csv_tunnel_profiles:
    :param snapshots:
    :param selected_tunnels:
    :param range:
    :return:
    """
    profiles = [line.split(sep=',') for line in open(csv_tunnel_profiles, 'r').readlines()[1:-1]]
    selected_dist = []
    selected_rads = []

    for snap, tunnel in zip(snapshots, selected_tunnels):
        for line in profiles:
            if line[0].strip() == snap and int(line[1]) == tunnel:
                if line[12].strip() == 'Z':
                    distances = [float(rad) for rad in line[13:]]
                    bottom = [i for i, dist in enumerate(distances) if dist < range[0]][0]
                    top = [i for i, dist in enumerate(distances) if dist > range[1]][0]
                    print('{} {}'.format(bottom, top))
                    # print(distances)
                    selected_dist.append(distances[bottom:top])
                if line[12].strip() == 'R':
                    rads = [float(rad) for rad in line[13:]]
                    rads = rads[bottom:top]
                    selected_rads.append(rads)
    return [selected_dist, selected_rads]


def plot_profile(parseds_profs, snapshots, outfile):
    """
    :param parseds_profs: list of parsed rdf profiles
    :return: plots parsed rdf profiles (no return)
    """
    mmplt.plot_simple_multiple(parseds_profs[1], parseds_profs[0], \
                               "pore radius [A]", "pore axis [A]", snapshots, outfile)


parser = argparse.ArgumentParser()
parser.add_argument("-tc", "--tunnel_chars", help="caver 3.0 tunnel characteristics in .csv format")
parser.add_argument("-tp", "--tunnel_profiles", help="caver 3.0 tunnel profiles in .csv format")
parser.add_argument("-s", "--pdbs", nargs='+', help="selected structures to plot")
parser.add_argument("-r", "--range", nargs='+', type=float, help="range of pore axis to plot")
parser.add_argument("-o", "--outfile", help="plot outfile name")
parser.add_argument("-p", "--profiles", nargs='+', type=int, help="selected profiles to plot")
parser.add_argument("-l", "--labels", nargs='+', help="optional labels")
args = parser.parse_args()


if args.profiles:
    print("info: parsing selected tunnels")
    parseds_profs = parse_selected_profiles(args.tunnel_profiles, args.pdbs, args.profiles, args.range)
    if args.labels:
        print("info: plotting using custom labels")
        labels = args.labels
    else:
        print("info: plotting using auto labels")
        labels = [snap + str(profile) for snap, profile in zip(args.pdbs, args.profiles)]

    plot_profile(parseds_profs, labels, args.outfile)

else:
    print("info: checking for reasonable tunnels")
    tunnel_chars = (parse_tunnel_characteristics(args.tunnel_chars))
    [pdbs, profiles] = find_main_tunnels(tunnel_chars)
    # print("info: parsing best tunnels")
    # parseds_profs = parse_selected_profiles(args.tunnel_profiles, pdbs, profiles, args.range)
    # print("info: plotting using standard labels")
    # labels = [snap + str(profile) for snap, profile in zip(pdbs, profiles)]

