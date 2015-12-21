# python 3
# plots vmd rmsd dat files
# michaladammichalowski@gmail.com
# 21.12.15 - creation
# EXAMPLE CALL: python3 mm_analysis_plotrmsd.py --rmsd_files rmsd_LONG.dat rmsd_STANDARD.dat --labels "prolonged equlibration" "standard conditions" --out_file rmsd_plot
# EXAMPLE CALL2: python3 mm_analysis_plotrmsd.py --rmsd_files rmsd_STANDARD_PRODUCTION.dat --labels  "standard conditions" --out_file rmsd_time_plot --timestep 2.0 --frequency 25000

import argparse
import numpy as np
import mm_lib_plots as mmplt


def read_rmsd(rmsd_files, timestep=False, frequency=False):
    """
    :param rmsd_files: rmsd files to read
    :param timestep: optional timestep for real time axis
    :param frequency: optional dcd frequency for real time axis
    :return: list of rmsd profiles
    """
    step_rmsd = []  # container for distances x
    vale_rmsd = []  # container for densities y

    for rmsd_file in rmsd_files:
        rmsd = open(rmsd_file).readlines()
        par_rmsd = [[float(col) for col in row.split()] for row in rmsd[2:-1]]  # first: create list rows

        if timestep:
            step_rmsd.append(np.arange(0, timestep*frequency*0.000001*len(par_rmsd), timestep*frequency*0.000001))
        else:
            step_rmsd.append([row[0] for row in par_rmsd])  # second: select x columns from rows and add to container
        vale_rmsd.append([row[1] for row in par_rmsd])  # second: select y columns from rows and add to container

    return [step_rmsd, vale_rmsd]


def plot_rdf(parseds_rmsd, labels, out_file, timestep=False):
    """
    :param parseds_rmsd: list of parsed rmsd profiles
    :param labels: data series labels
    :param out_file: file to save plot
    :param timestep: send timestep for real time label
    :return: plots parsed rdf profiles (no return)
    """
    if timestep:
        x_label = "time [ns]"
    else:
        x_label = "step []"
    mmplt.plot_simple_multiple(parseds_rmsd[0], parseds_rmsd[1], x_label, "RMSD [A]", labels, out_file)

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--rmsd_files", nargs='+', help="rmsd mfiles to plot")
parser.add_argument("-l", "--labels", nargs='+', help="labels for data series")
parser.add_argument("-o", "--out_file", help="outfile name")
parser.add_argument("-t", "--timestep", type=float, help="timestep")
parser.add_argument("-f", "--frequency", type=int, help="dcd save frequency")
args = parser.parse_args()

rmsds = read_rmsd(args.rmsd_files, args.timestep, args.frequency)
plot_rdf(rmsds, args.labels, args.out_file, args.timestep)
