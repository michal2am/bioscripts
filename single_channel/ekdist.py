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

def fit_ekdist_openings(scns, tres, taus, areas, color, plot_name, lim):

    rec = ekrecord.SingleChannelRecord()
    rec.load_SCN_file(scns)
    rec.tres = tres
    print(rec)
    #ekplot.plot_xlog_interval_histogram(rec.opint, rec.tres, shut=False)
    #plt.show()


    expPDF = dceq.MultiExponentialPDF(np.asarray(rec.opint), taus=np.asarray(taus), areas=np.asarray(areas))
    theta = expPDF.theta
    print('Start LogLikelihood =', expPDF.loglik(theta))
    res = minimize(expPDF.loglik, theta, method='Nelder-Mead')
    print(res)
    expPDF.theta = res.x

    ekplot.plot_xlog_interval_histogram_fit(rec.opint, rec.tres, expPDF.to_plot, res.x, shut=False)
    print(expPDF)

    ax = plt.gca()
    for line in ax.get_lines():
        line.set_color(color)

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(12)

    for axis in ['top', 'right']:
        ax.spines[axis].set_linewidth(0)
    ax.tick_params(width=2, which='both')
    ax.set_xlim([3e-5, lim])

    plt.tight_layout()
    plt.savefig(plot_name, dpi=600)
    plt.show()

def fit_ekdist_shuts(scns, tres, taus, areas, color, plot_name, lim):

    rec = ekrecord.SingleChannelRecord()
    rec.load_SCN_file(scns)
    rec.tres = tres
    print(rec)
    # ekplot.pl/ot_xlog_interval_histogram(rec.shint, rec.tres, shut=True)
    # plt.show()

    expPDF = dceq.MultiExponentialPDF(np.asarray(rec.shint), taus=np.asarray(taus), areas=np.asarray(areas))
    theta = expPDF.theta
    print('Start LogLikelihood =', expPDF.loglik(theta))
    res = minimize(expPDF.loglik, theta, method='Nelder-Mead')
    print(res)
    expPDF.theta = res.x

    ekplot.plot_xlog_interval_histogram_fit(rec.shint, rec.tres, expPDF.to_plot, res.x, shut=True)
    print(expPDF)
    ax = plt.gca()
    for line in ax.get_lines():
        line.set_color(color)

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(12)

    for axis in ['top', 'right']:
        ax.spines[axis].set_linewidth(0)
    ax.tick_params(width=2, which='both')
    ax.set_xlim([3e-5, lim])

    plt.tight_layout()
    plt.savefig(plot_name, dpi=600)
    plt.show()


fit_ekdist_shuts(['53WT3_K1.SCN'], 40e-6, [0.02e-3, 0.156e-3, 0.55e-3, 6.04e-3], [0.67, 0.25, 0.07, 0.006], '#0070c0', 'V53WT3_shuts.png', 1e-1)
fit_ekdist_openings(['53WT3_K1.SCN'], 40e-6, [1.28e-3, 2.19e-3], [0.48, 0.52], '#ff0000', 'V53WT3_openings.png', 1e-1)



#fit_ekdist_openings(scns, tres, taus_op, areas_op)
#fit_ekdist_shuts(scns, tres, taus_sh, areas_sh)

#fit_ekdist_shuts(['V53A2_K2.SCN', 'V53A2_K4.SCN', 'V53A2_K8.SCN'], 50e-6, [0.4e-3, 2.44e-3, 33.67e-3], [0.65, 0.33, 0.01])
#fit_ekdist_shuts(['V53A2_K2.SCN', 'V53A2_K4.SCN', 'V53A2_K8.SCN'], 50e-6, [0.04e-3, 0.26e-3, 1.56e-3, 21.72e-3], [0.59, 0.31, 0.09, 0.01])
#fit_ekdist_shuts(['V53A2_K2.SCN', 'V53A2_K4.SCN', 'V53A2_K8.SCN'], 50e-6, [0.0005, 0.002, 0.001, 0.01], [0.25, 0.25, 0.25, 0.25])

#fit_ekdist_shuts(['V53A3_K6.SCN'], 50e-6, [0.53e-3, 1.49e-3, 20.12e-3], [0.78, 0.21,0.01])
#fit_ekdist_shuts(['V53A3_K6.SCN'], 50e-6, [0.0005, 0.002, 0.001, 0.01], [0.25, 0.25, 0.25, 0.25])

#fit_ekdist_shuts(['V53A8_K2.SCN', 'V53A8_K3.SCN', 'V53A8_K6.SCN'], 50e-6, [0.64e-3, 2.74e-3, 25.91e-3], [0.74, 0.23, 0.02])
#fit_ekdist_shuts(['V53A8_K2.SCN', 'V53A8_K3.SCN', 'V53A8_K6.SCN'], 50e-6, [0.0005, 0.002, 0.001, 0.01], [0.25, 0.25, 0.25, 0.25])

#fit_ekdist_shuts(['V53A9K14.SCN', 'V53A9K17.SCN'], 50e-6, [0.60e-3, 2.32e-3, 35.12e-3], [0.69, 0.30, 0.01])
#fit_ekdist_shuts(['V53A9K14.SCN', 'V53A9K17.SCN'], 50e-6, [0.0005, 0.002, 0.001, 0.01], [0.25, 0.25, 0.25, 0.25])

# paper plots:

#fit_ekdist_shuts(['V53A2_K2.SCN', 'V53A2_K4.SCN', 'V53A2_K8.SCN'], 50e-6, [0.4e-3, 2.44e-3, 33.67e-3], [0.65, 0.33, 0.01], '#dd8452', 'V53A_shuts.png', 1e0)
#fit_ekdist_openings(['V53A2_K2.SCN', 'V53A2_K4.SCN', 'V53A2_K8.SCN'], 50e-6, [0.58e-3, 1.92e-3], [0.88, 0.12], '#dd8452', 'V53A_openings.png', 1e-1)

#fit_ekdist_shuts(['53E4_K1.SCN', '53E4_K2.SCN', '53E4_K3.SCN'], 40e-6, [2.64e-3, 11.76e-3], [0.66, 0.34], '#55a868', 'V53E_shuts.png', 1e0)
#fit_ekdist_openings(['53E4_K1.SCN', '53E4_K2.SCN', '53E4_K3.SCN'], 40e-6, [0.18e-3, 0.86e-3], [0.89, 0.11], '#55a868', 'V53E_openings.png', 1e-2)

#fit_ekdist_shuts(['53H3_K29.SCN', '53H3_K42.SCN', '53H3_K58.SCN', '53H3_K39.SCN'], 40e-6, [0.26e-3, 1.58e-3, 14.11e3], [0.61, 0.34, 0.05], '#4c72b0', 'V53H_shuts.png', 1e0)
#fit_ekdist_openings(['53H3_K29.SCN', '53H3_K42.SCN', '53H3_K58.SCN', '53H3_K39.SCN'], 40e-6, [0.56e-3, 1.77e-3], [0.81, 0.19], '#4c72b0', 'V53H_openings.png', 1e-1)

#fit_ekdist_shuts(['53K3_K1.SCN'], 50e-6, [0.04e-3, 0.56e-3, 2.46e-3, 26.19e-3], [0.60, 0.20, 0.16, 0.04], '#c44e52', 'V53K_shuts.png', 1e0)
#fit_ekdist_shuts(['53K3_K1.SCN'], 50e-6, [0.042e-3, 0.56e-3, 2.46e3, 26.19e3], [0.60, 0.20, 0.16, 0.04], '#c44e52', 'V53K_shuts.png', 1e0)


#basic scenario
#fit_ekdist_shuts(['53WT2_K1.SCN','53WT2_K2.SCN','53WT2_K3.SCN','53WT2_K4.SCN'], 50e-6, [0.02e-3, 0.2e-3, 1.08e-3, 8.84e-3],[0.64, 0.28, 0.07, 0.007], '#8c8c8c', 'V53WT_shuts.png', 1e0)
#fit_ekdist_openings(['53WT2_K1.SCN','53WT2_K2.SCN','53WT2_K3.SCN','53WT2_K4.SCN'], 50e-6, [1.28e-3, 3.20e-3], [0.55, 0.45], '#8c8c8c', 'V53WT_openings.png', 1e-1)
#WT3
#fit_ekdist_shuts(['53WT3_K1.SCN'], 50e-6, [0.02e-3, 0.15e-3, 0.54e-3, 6.05e-3],[0.67, 0.25, 0.07, 0.007], '#8c8c8c', 'V53WT_shuts.png', 1e0)


#A dd8452
#K 55a868
#H 4c72b0
#K c44e52
#WT 8c8c8c