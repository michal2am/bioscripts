import pandas as pd
from dcpyps import dcio
import seaborn as sns
import matplotlib.pyplot as plt
import pyabf
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--scn_file')

args = parser.parse_args()

# SCN stuff:

scn_file_name = args.scn_file

header = dcio.scn_read_header(scn_file_name)
scan = dcio.scn_read_data(scn_file_name, header[3])     # header[3] seems to contain needed dict with the parameters
scan = pd.DataFrame({'period_scan': scan[0], 'scan_amp': scan[1], 'scan_cat': scan[-1]})
print(scan)

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