import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from dcpyps import dcio
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--clmpf_files', nargs='+')
args = parser.parse_args()


for clmpf in args.clmpf_files:
    clmpf_dwells = pd.read_excel(clmpf, header=None, names=['Trace', 'Serach', 'Level', 'State', 'Start', 'End', 'Amplitude', 'Amplitude_SD', 'Dwell_time', 'Peak_to_peak', 'Interevent'])
    #print()
    #print(clmpf_dwells)
    #dcio.scn_write(clmpf_dwells.Dwell_time.values.tolist(), clmpf_dwells.Level.values.tolist(), len(clmpf_dwells.Dwell_time.values.tolist())*[0], filename=clmpf.split('.')[0] + '_clmpf.SCN')
    print(clmpf_dwells.Dwell_time.values.tolist())
    print(clmpf_dwells.Amplitude.values.tolist())
    print(len(clmpf_dwells.Dwell_time.values.tolist())*[0])
    # why this stopped to work???? probably did not work ever before
    dcio.scn_write(clmpf_dwells.Dwell_time.values.tolist(), clmpf_dwells.Amplitude.values.tolist(), len(clmpf_dwells.Dwell_time.values.tolist())*[0], filename=clmpf.split('.')[0] + '_cfA.SCN')


    #header = dcio.scn_read_header(clmpf.split('.')[0] + '_clmpf.SCN')
    #print(header)

#clampfit_dwell_11 = pd.read_csv('WT_25_C1_Clampfit_11.csv')
#g = sns.relplot(kind='line', x='Start' , y='Level', data=clampfit_dwell_11, drawstyle='steps-pre', color='grey')
#plt.show()



#print(clampfit_dwell_11)

#dcio.scn_write(clampfit_dwell_11.Dwell_time.values.tolist(), clampfit_dwell_11.Level.values.tolist(), len(clampfit_dwell_11.Dwell_time.values.tolist())*[0], filename='WT_25_C1_Clampfit_11.SCN')
#header = dcio.scn_read_header('WT_25_C1_Clampfit_11.SCN')
#print(header)
#scan = dcio.scn_read_data('WT_25_C1_Clampfit_11.SCN', header[3])
#print(scan)
#ioffset, nint, calfac2, header = dcio.scn_read_header('f45.SCN')
#print(type(header['ioffset']))
#dcio.scn_read_data('f45.SCN', dcio.scn_read_header('f45.SCN'))

