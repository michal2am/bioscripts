# python 3
# plots vmd rmsd dat files
# michaladammichalowski@gmail.com
# 21.12.15 - creation
# 10.03.15 - refactor
# EXAMPLE CALL: python3 mm_analysis_plotrmsd.py --rmsd_files rmsd_LONG.dat rmsd_STANDARD.dat --labels "prolonged equlibration" "standard conditions" --out_file rmsd_plot --stabper 500
# EXAMPLE CALL2: python3 mm_analysis_plotrmsd.py --rmsd_files rmsd_STANDARD_PRODUCTION.dat --labels  "standard conditions" --out_file rmsd_time_plot --timestep 2.0 --frequency 25000 --stabper 500

import argparse
import numpy as np
from scipy import stats
import mm_lib_plots as mmplt


class RMSD:

    def __init__(self, in_files, out_file, labels, timestep=False, frequency=False, period=False, speriod=False):
        """
        :param in_files: rmsd files to read
        :param out_file: file to save plot
        :param timestep: optional timestep for real time axis
        :param frequency: optional dcd frequency for real time axis
        :param period:
        :param labels: data series labels
        :return:
        """
        self.in_files = in_files
        self.out_file = out_file
        self.timestep = timestep
        self.frequency = frequency
        self.period = period
        self.speriod = speriod
        self.labels = labels
        self.rmsds = self.read_rmsd()

    def read_rmsd(self):
        """
        :return: list of rmsd profiles
        """
        rmsds = list()

        for rmsd_file in self.in_files:

            par_rmsd = np.genfromtxt(rmsd_file, skip_header=2)

            if self.timestep:
                par_rmsd[:, 0] *= self.timestep*self.frequency*0.000001

            if self.period:
                par_rmsd = par_rmsd[self.period:]

            rmsds.append(par_rmsd)
        return rmsds

    def check_stab(self):
        """
        :return: calculates linear regression for selected period and first derivative for whole trajectory
        """

        for i, rmsd in enumerate(self.rmsds):

            slope, intercept, r_value, p_value, std_err = stats.linregress(rmsd[self.speriod:, 0], rmsd[self.speriod:, 1])
            print("DATA: {} Slope: {:.4f} Intercept: {:.4f} R value: {:.3f} P value: {:.3f}"
                  .format(self.labels[i], slope, intercept, r_value, p_value))
            fit = rmsd[:, 0:1] * slope + intercept  # select vs slice notation here to get column vector
            self.rmsds[i] = np.append(rmsd, fit, axis=1)  # enumeration for list element replacement

            deriv = np.append(np.zeros(1), np.diff(rmsd[:, 1])) * 3  # magic scaling factor
            deriv.shape = (len(deriv), 1)
            self.rmsds[i] = np.append(self.rmsds[i], deriv, axis=1)

    def plot_rmsd(self):
        """
        :return: plots parsed rdf profiles
        """
        x_label = "time [ns]" if self.timestep else "step []"
        mmplt.plot_simple_multiple_numpy(self.rmsds, x_label, "RMSD [A]", self.labels, self.out_file, sizex=6, sizey=6)

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--rmsd_files", nargs='+', help="rmsd files to plot")
parser.add_argument("-l", "--labels", nargs='+', help="labels for data series")
parser.add_argument("-o", "--out_file", help="outfile name")
parser.add_argument("-t", "--timestep", type=float, help="timestep")
parser.add_argument("-f", "--frequency", type=int, help="dcd save frequency")
parser.add_argument("-p", "--period", type=int, help="custom start frame")
parser.add_argument("-s", "--stabper", type=float, help="stability check start")
args = parser.parse_args()

rmsds = RMSD(args.rmsd_files, args.out_file, args.labels, args.timestep, args.frequency, args.period, args.stabper)
rmsds.check_stab()
rmsds.plot_rmsd()
