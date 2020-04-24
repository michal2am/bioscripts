import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from dcpyps import dcio

input = '89145k1.csv'

sch = pd.read_csv(input, header=None, usecols=[0, 1], names=['time_real', 'current_real'])
sch.loc[:, 'current_real'] = sch.loc[:, 'current_real']/10 + 2.4

header = dcio.scn_read_header('89145K1.SCN')
#print(header[3]['ioffest'])

scan = dcio.scn_read_data('89145K1.SCN', dcio.scn_read_header('89145K1.SCN')[3])
scan = pd.DataFrame({'period_scan':scan[0], 'scan_amp':scan[1], 'scan_cat': scan[-1]})

scan.loc[:, 'time_scan'] = scan.loc[:, 'period_scan'].cumsum()
scan.loc[:, 'time_scan'] += 0.0895

scan.loc[:, 'scan_amp_n'] = ((scan.loc[:, 'scan_amp'] - min(scan.loc[:, 'scan_amp']))/min(scan.loc[:, 'scan_amp'])) +2.1

scan.loc[:, 'scan_cat_w'] = np.NaN
scan.loc[:, 'scan_cat_w'][scan.loc[:, 'scan_cat'] > 0] = 2


# scan.loc[:, 'scan_amp_n'][scan.loc[:, 'scan_amp'] < -1] = 2.1
# scan.loc[:, 'scan_amp_n'][scan.loc[:, 'scan_amp'] > -0.5] = 1.1



output = 'result_89145k1_2kFilter_Red.csv'


fit = pd.read_csv(output, usecols=['Time', 'prediction'])
fit = pd.DataFrame(fit)
fit.rename(columns={'Time': 'time_deep'}, inplace=True)
fit.loc[:, 'prediction'] = fit.loc[:, 'prediction'] - 2
fit.loc[:, 'time_deep'] -= 0.0328


'''
output2 = 'result_89145k1_2kFilter.csv'

fit2 = pd.read_csv(output2, usecols=['Time', 'prediction'])
fit2 = pd.DataFrame(fit2)
fit2.rename(columns={'Time': 'time_deep_2', 'prediction': 'prediction_2'}, inplace=True)
print(fit2)
print(fit2.loc[:, 'prediction_2'])
fit2.loc[:, 'prediction_2'] = fit2.loc[:, 'prediction_2'] - 3.1
print('dupa2')
fit2.loc[:, 'time_deep_2'] -= 0.033
print(fit2)
'''

sch_fit = pd.concat([sch, fit, scan], axis=1)

print(sch_fit)
sch_fit.to_csv('results.csv')

# sch_fit = pd.concat([sch, fit, fit2], axis=1)
g = sns.relplot(kind='line', x="time_real", y="current_real", data=sch_fit, color='grey')
g.map(sns.lineplot, x='time_scan', y='scan_amp_n', data=scan, color='teal', drawstyle='steps-pre')
# g.map(sns.lineplot, x='time_scan', y='scan_amp_n', data=scan[scan.loc[:, 'scan_cat'] < 1.9], color='turquoise', drawstyle='steps-pre')
g.map(sns.scatterplot, x='time_scan', y='scan_cat_w', data=scan, color='tomato')
g.map(sns.lineplot, x='time_deep', y='prediction', data=sch_fit, color='goldenrod', drawstyle='steps-pre')

# g.map(sns.lineplot, x='time_deep_2', y='prediction_2', data=sch_fit, color='pink')

plt.show()
