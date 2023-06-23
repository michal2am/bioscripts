# python ~/repos/bioscripts/single_channel/ideal\&files/scn2csv.py --scn_file WT_25_C1.SCN --out_file WT_25_C1.csv --abf_file WT_25_c1.abf --plot_abf_scn yes

import pandas as pd
from dcpyps import dcio
import seaborn as sns
import matplotlib.pyplot as plt
import pyabf
import  argparse

parser = argparse.ArgumentParser()
parser.add_argument('--scn_file')
parser.add_argument('--out_file')
# those are optional, only for scn. vs abf plot
parser.add_argument('--abf_file')
parser.add_argument('--plot_abf_scn')

args = parser.parse_args()

scn_file_name = args.scn_file

header = dcio.scn_read_header(scn_file_name)
scan = dcio.scn_read_data(scn_file_name, header[3])     # header[3] seems to contain needed dict with the parameters
scan = pd.DataFrame({'period_scan': scan[0], 'scan_amp': scan[1], 'scan_cat': scan[-1]})


"""    
    0 = all OK;
    1 = amplitude dubious = bit 0;
    2 = amplitude fixed = bit 1;
    4 = amplitude of opening constrained (see fixamp) = bit 2;
    8 = duration unusable = bit 3; etc
"""

# print(scan[scan['scan_cat'] == 8])

# normalized amplitude
# scan.loc[:, 'scan_amp_n'] = ((scan.loc[:, 'scan_amp'] - min(scan.loc[:, 'scan_amp']))/min(scan.loc[:, 'scan_amp']))
# just binary 0 or 1 amplitudes
scan.loc[scan['scan_amp'] != 0, 'scan_amp_b'] = 1
scan.loc[scan['scan_amp'] == 0, 'scan_amp_b'] = 0
scan['period_scan'] = scan['period_scan']*1000

# adding up series of neighbouring shut or open dwell times
start_event_no = 0
scan_shifted = scan['scan_amp_b'].shift()
scan['scan_amp_b_prev'] = scan_shifted


def event_no(row):
    global start_event_no
    if row['scan_amp_b'] != row['scan_amp_b_prev']:
        start_event_no += 1
        return start_event_no
    else:
        return start_event_no


scan['event_no'] = scan.apply(lambda x: event_no(x), axis=1)

scan_dwell = scan.groupby('event_no').agg({'period_scan': 'sum', 'scan_amp_b': 'max'})
scan_dwell['time_scan'] = scan_dwell['period_scan'].cumsum()
scan_dwell.loc[scan_dwell['period_scan'] < 0.05, 'under_res'] = -0.1
scan_dwell.loc[scan_dwell['period_scan'] >= 0.05, 'under_res'] = -0.2
print(scan_dwell)
scan_dwell.to_csv(args.out_file, index=False)

clampfit_dwell_09 = pd.read_csv('WT_25_C1_Clampfit_09.csv')
clampfit_dwell_07 = pd.read_csv('WT_25_C1_Clampfit_07.csv')
clampfit_dwell_05 = pd.read_csv('WT_25_C1_Clampfit_05.csv')

clampfit_dwell_095 = pd.read_csv('WT_25_C1_Clampfit_095.csv')
clampfit_dwell_10 = pd.read_csv('WT_25_C1_Clampfit_10.csv')
clampfit_dwell_105 = pd.read_csv('WT_25_C1_Clampfit_105.csv')
clampfit_dwell_11 = pd.read_csv('WT_25_C1_Clampfit_11.csv')

clampfit_dwell_10_20 = pd.read_csv('WT_25_C1_Clampfit_10_20p.csv')
clampfit_dwell_11_20 = pd.read_csv('WT_25_C1_Clampfit_11_20p.csv')

# TODO: 0.95 == 1.6492

#print(clampfit_dwell)

if args.plot_abf_scn:

    # upward currents are openings

    abf_time_shift = 0.05965
    abf_amplitude = 2.5
    abf_amplitude_shift = 0

    abf = pyabf.ABF(args.abf_file)
    abf.setSweep(0)
    #plt.plot(abf.sweepX, abf.sweepY)
    #plt.show()

    # below some magic numbers to scale amplitude and shift time, this should be read from scn/prt
    g = sns.relplot(kind='line', x="time_scan", y="scan_amp_b", data=scan_dwell, drawstyle='steps-pre', color='grey')
    g.map(sns.lineplot, x=scan_dwell.time_scan, y=scan_dwell.under_res, drawstyle='steps-pre', color='grey')
    # abf recording:
    g.map(sns.lineplot, x=(abf.sweepX - abf_time_shift)*1000, y=((abf.sweepY - abf_amplitude_shift )/abf_amplitude))
    # clampfit ideal:
    # g.map(sns.lineplot, x=clampfit_dwell_05.Start - 59.7, y=(clampfit_dwell_05.Level/10) + 1.05, drawstyle='steps-post', color='red')
    # g.map(sns.lineplot, x=clampfit_dwell_07.Start - 59.7, y=(clampfit_dwell_07.Level/10) + 1.16, drawstyle='steps-post', color='red')
    # g.map(sns.lineplot, x=clampfit_dwell_09.Start - 59.7, y=(clampfit_dwell_09.Level/10) + 1.27, drawstyle='steps-post', color='red')

    #g.map(sns.lineplot, x=clampfit_dwell_095.Start - 59.7, y=(clampfit_dwell_095.Level/10) + 1.38, drawstyle='steps-post', color='red')
    g.map(sns.lineplot, x=clampfit_dwell_10.Start - 59.7, y=(clampfit_dwell_10.Level/10) + 1.49, drawstyle='steps-post', color='red')
    g.map(sns.lineplot, x=clampfit_dwell_10_20.Start - 59.7, y=(clampfit_dwell_10_20.Level / 10) + 1.60,drawstyle='steps-post', color='red')
    #g.map(sns.lineplot, x=clampfit_dwell_105.Start - 59.7, y=(clampfit_dwell_105.Level/10) + 1.60, drawstyle='steps-post', color='red')
    g.map(sns.lineplot, x=clampfit_dwell_11.Start - 59.7, y=(clampfit_dwell_11.Level/10) + 1.71, drawstyle='steps-post', color='red')
    g.map(sns.lineplot, x=clampfit_dwell_11_20.Start - 59.7, y=(clampfit_dwell_11_20.Level/10) + 1.82, drawstyle='steps-post', color='red')


    plt.show()