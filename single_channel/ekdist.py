import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from dcpyps.ekdist import ekrecord
from dcpyps.ekdist import ekplot
from dcpyps import dcequations as dceq
from dcpyps import dataset
import numpy as np
from scipy.optimize import minimize, bisect


rec = ekrecord.SingleChannelRecord()
#TODO: should work for multiple scns, copy from SCALCS dataset.SCRecord
rec.load_SCN_file(['V53A2_K2.SCN', 'V53A2_K4.SCN', 'V53A2_K8.SCN'])
rec.tres = 50e-6
print(rec)
#ekplot.plot_fitted_amplitude_histogram(rec, fc, n=2)
ekplot.plot_xlog_interval_histogram(rec.opint, rec.tres, shut=False)
plt.show()
ekplot.plot_xlog_interval_histogram(rec.shint, rec.tres, shut=True)
plt.show()


taus = [0.0005, 0.002]
areas = [0.5, 0.5]
expPDF = dceq.MultiExponentialPDF(np.asarray(rec.opint),
                                         taus=np.asarray(taus), areas=np.asarray(areas))
theta = expPDF.theta
print('Start LogLikelihood =', expPDF.loglik(theta))
res = minimize(expPDF.loglik, theta, method='Nelder-Mead')
print(res)
expPDF.theta = res.x

ekplot.plot_xlog_interval_histogram_fit(rec.opint, rec.tres, expPDF.to_plot, res.x, shut=False)
print(expPDF)
plt.show()