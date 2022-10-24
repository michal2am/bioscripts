import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from dcpyps.ekdist import ekrecord
from dcpyps.ekdist import ekplot
from dcpyps import dcequations as dceq
from dcpyps import dataset
import numpy as np
from scipy.optimize import minimize, bisect


taus_op = [0.0005, 0.002]
areas_op = [0.5, 0.5]
taus_sh = [0.0005, 0.002, 0.001, 0.01]
areas_sh = [0.25, 0.25, 0.25, 0.25]
scns = ['V53A2_K2.SCN', 'V53A2_K4.SCN', 'V53A2_K8.SCN']
tres = 50e-6

def fit_ekdist_openings(scns, tres, taus, areas):

    rec = ekrecord.SingleChannelRecord()
    rec.load_SCN_file(scns)
    rec.tres = tres
    print(rec)
    ekplot.plot_xlog_interval_histogram(rec.opint, rec.tres, shut=False)
    plt.show()


    expPDF = dceq.MultiExponentialPDF(np.asarray(rec.opint), taus=np.asarray(taus), areas=np.asarray(areas))
    theta = expPDF.theta
    print('Start LogLikelihood =', expPDF.loglik(theta))
    res = minimize(expPDF.loglik, theta, method='Nelder-Mead')
    print(res)
    expPDF.theta = res.x

    ekplot.plot_xlog_interval_histogram_fit(rec.opint, rec.tres, expPDF.to_plot, res.x, shut=False)
    print(expPDF)
    plt.show()

def fit_ekdist_shuts(scns, tres, taus, areas):

    rec = ekrecord.SingleChannelRecord()
    rec.load_SCN_file(scns)
    rec.tres = tres
    print(rec)
    # ekplot.plot_xlog_interval_histogram(rec.shint, rec.tres, shut=True)
    # plt.show()

    expPDF = dceq.MultiExponentialPDF(np.asarray(rec.shint), taus=np.asarray(taus), areas=np.asarray(areas))
    theta = expPDF.theta
    print('Start LogLikelihood =', expPDF.loglik(theta))
    res = minimize(expPDF.loglik, theta, method='Nelder-Mead')
    print(res)
    expPDF.theta = res.x

    ekplot.plot_xlog_interval_histogram_fit(rec.shint, rec.tres, expPDF.to_plot, res.x, shut=True)
    print(expPDF)
    plt.show()

#fit_ekdist_openings(scns, tres, taus_op, areas_op)
#fit_ekdist_shuts(scns, tres, taus_sh, areas_sh)

fit_ekdist_shuts(['V53A2_K2.SCN', 'V53A2_K4.SCN', 'V53A2_K8.SCN'], 50e-6, [0.4e-3, 2.44e-3, 33.67e-3], [0.65, 0.33, 0.01])
#fit_ekdist_shuts(['V53A2_K2.SCN', 'V53A2_K4.SCN', 'V53A2_K8.SCN'], 50e-6, [0.04e-3, 0.26e-3, 1.56e-3, 21.72e-3], [0.59, 0.31, 0.09, 0.01])
fit_ekdist_shuts(['V53A2_K2.SCN', 'V53A2_K4.SCN', 'V53A2_K8.SCN'], 50e-6, [0.0005, 0.002, 0.001, 0.01], [0.25, 0.25, 0.25, 0.25])

fit_ekdist_shuts(['V53A3_K6.SCN'], 50e-6, [0.53e-3, 1.49e-3, 20.12e-3], [0.78, 0.21,0.01])
fit_ekdist_shuts(['V53A3_K6.SCN'], 50e-6, [0.0005, 0.002, 0.001, 0.01], [0.25, 0.25, 0.25, 0.25])

fit_ekdist_shuts(['V53A8_K2.SCN', 'V53A8_K3.SCN', 'V53A8_K6.SCN'], 50e-6, [0.64e-3, 2.74e-3, 25.91e-3], [0.74, 0.23, 0.02])
fit_ekdist_shuts(['V53A8_K2.SCN', 'V53A8_K3.SCN', 'V53A8_K6.SCN'], 50e-6, [0.0005, 0.002, 0.001, 0.01], [0.25, 0.25, 0.25, 0.25])

fit_ekdist_shuts(['V53A9K14.SCN', 'V53A9K17.SCN'], 50e-6, [0.60e-3, 2.32e-3, 35.12e-3], [0.69, 0.30, 0.01])
fit_ekdist_shuts(['V53A9K14.SCN', 'V53A9K17.SCN'], 50e-6, [0.0005, 0.002, 0.001, 0.01], [0.25, 0.25, 0.25, 0.25])