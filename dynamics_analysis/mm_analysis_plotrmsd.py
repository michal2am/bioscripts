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
        :param period: cut data from
        :param speriod: approximate data from
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

        if self.timestep: self.timeframe = self.timestep*self.frequency*0.000001
        self.rmsds = self.read_rmsd()[0]
        self.approx_sta = self.read_rmsd()[1]
        self.approx_len = self.read_rmsd()[2]

    def read_rmsd(self):
        """
        :return: list of rmsd profiles
        """
        rmsds = list()

        for rmsd_file in self.in_files:

            par_rmsd = np.genfromtxt(rmsd_file, skip_header=2)

            if self.period:
                par_rmsd = par_rmsd[self.period[0]:self.period[1]]

            if self.timestep:
                par_rmsd[:, 0] *= self.timeframe

            rmsds.append(par_rmsd)

            approx_sta = self.speriod*self.timeframe if self.timestep else self.speriod
            approx_len = rmsds[0][-1][0] - approx_sta

        return rmsds, approx_sta, approx_len

    def check_stab(self):
        """
        :return: calculates linear regression for selected period and first derivative for whole trajectory
        """

        if self.timestep:
            print('Stability check from step {:.0f} equal to {:.2f} ns'.format(self.speriod, self.approx_sta))
            print('Stability period length {:.2f} ns'.format(self.approx_len))


        print("fit                  slope    intercept R-value P-value")
        print("                     [A^2/ns] [A^2]")


        for i, rmsd in enumerate(self.rmsds):

            slope, intercept, r_value, p_value, std_err = stats.linregress(rmsd[self.speriod:, 0], rmsd[self.speriod:, 1])
            print("{: <20} {:.5f}  {:.2f}      {:.3f}   {:.3f}"
                  .format(self.labels[i], slope, intercept, r_value, 1 - p_value))
            fit = rmsd[:, 0:1] * slope + intercept  # select vs slice notation here to get column vector
            self.rmsds[i] = np.append(rmsd, fit, axis=1)  # enumeration for list element replacement

            if args.deriv:
                deriv = np.append(np.zeros(1), np.diff(rmsd[:, 1])) * 3  # magic scaling factor
                deriv.shape = (len(deriv), 1)
                self.rmsds[i] = np.append(self.rmsds[i], deriv, axis=1)

    def plot_rmsd(self):
        """
        :return: plots parsed rdf profiles
        """
        x_label = "time [ns]" if self.timestep else "step []"
        mmplt.plot_simple_multiple_numpy(self.rmsds, x_label, "RMSD [A]", self.labels, self.out_file,
                                         sizex=1.5, sizey=1.5, ylimit=[args.ranges[0], args.ranges[1]],
                                         back=[[(self.approx_sta, self.approx_len)],
                                               (args.ranges[0], args.ranges[1] - args.ranges[0])])

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--rmsd_files", nargs='+', help="rmsd files to plot")
parser.add_argument("-l", "--labels", nargs='+', help="labels for data series")
parser.add_argument("-o", "--out_file", help="outfile name")
parser.add_argument("-t", "--timestep", type=float, help="timestep")
parser.add_argument("-f", "--frequency", type=int, help="dcd save frequency")
parser.add_argument("-p", "--period", nargs='+', type=int, help="custom start/stop frame")
parser.add_argument("-s", "--stabper", type=float, help="stability check start frame")
parser.add_argument("-d", "--deriv", help="add derivative")
parser.add_argument("-a", "--ranges", nargs='+', type=float, help="RMSD plot range")
args = parser.parse_args()

rmsds = RMSD(args.rmsd_files, args.out_file, args.labels, args.timestep, args.frequency, args.period, args.stabper)
rmsds.check_stab()
rmsds.plot_rmsd()
