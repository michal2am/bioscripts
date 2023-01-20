import pandas as pd
import numpy as np
from dcpyps import dcio

scn_file_name = '89145K1.SCN'

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

# this transfers dwell times to absolute time,
#scan.loc[:, 'time_scan'] = scan.loc[:, 'period_scan'].cumsum()
# ???
#scan.loc[:, 'time_scan'] += 0.0895


scan.loc[:, 'scan_amp_n'] = ((scan.loc[:, 'scan_amp'] - min(scan.loc[:, 'scan_amp']))/min(scan.loc[:, 'scan_amp'])) +2.1

#scan.loc[:, 'scan_cat_w'] = np.NaN
#scan.loc[:, 'scan_cat_w'][scan.loc[:, 'scan_cat'] > 0] = 2

print(scan)