import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from dcpyps import dcio


clampfit_dwell_11 = pd.read_csv('WT_25_C1_Clampfit_11.csv')
g = sns.relplot(kind='line', x='Start' , y='Level', data=clampfit_dwell_11, drawstyle='steps-pre', color='grey')
#plt.show()



print(clampfit_dwell_11)

dcio.scn_write(clampfit_dwell_11.Dwell_time.values.tolist(), clampfit_dwell_11.Level.values.tolist(), len(clampfit_dwell_11.Dwell_time.values.tolist())*[0], filename='WT_25_C1_Clampfit_11.SCN')
header = dcio.scn_read_header('WT_25_C1_Clampfit_11.SCN')
print(header)
scan = dcio.scn_read_data('WT_25_C1_Clampfit_11.SCN', header[3])
print(scan)
#ioffset, nint, calfac2, header = dcio.scn_read_header('f45.SCN')
#print(type(header['ioffset']))
#dcio.scn_read_data('f45.SCN', dcio.scn_read_header('f45.SCN'))

