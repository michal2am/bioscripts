import pandas as pd
import numpy as np
from dcpyps import dcio
import seaborn as sns
import matplotlib.pyplot as plt
import pyabf

abf = pyabf.ABF("WT_25_c1.abf")
abf.setSweep(0)
#plt.plot(abf.sweepX, abf.sweepY)
#plt.show()


scn_file_name = 'WT_25_C1.SCN'

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

print(scan_dwell)
# below some magic numbers to scale amplitude and shift time, this should be read from scn/prt
g = sns.relplot(kind='line', x="time_scan", y="scan_amp_b", data=scan_dwell, drawstyle='steps-pre', color='grey')
g.map(sns.lineplot, x=(abf.sweepX - 0.0625)*1000, y=(abf.sweepY -0.0)/2.5)
plt.show()
