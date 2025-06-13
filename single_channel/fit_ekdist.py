import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from dcpyps.ekdist import ekrecord
from dcpyps.ekdist import ekplot
from dcpyps import dcequations as dceq
from dcpyps import dataset
import numpy as np
from scipy.optimize import minimize, bisect

def fit_from_config_4sh1op(ekdist_config):
    config = pd.read_csv(ekdist_config)
    results = []


    for cell in config.cell_id.unique():

        single_cell = config[config.loc[:, 'cell_id'] == cell].copy()
        single_cell.reset_index(inplace=True)

        # ideal time distribution areas from hjcfit works better for shut dwell time fitting
        # taus are rather similar in both cases

        sh_p1id = single_cell.at[0, 'sh_p1id']/100
        sh_p2id = single_cell.at[0, 'sh_p2id']/100
        sh_p3id = single_cell.at[0, 'sh_p3id']/100
        sh_p4id = single_cell.at[0, 'sh_p4id']/100

        sh_t1id = single_cell.at[0, 'sh_t1id']*1e-3
        sh_t2id = single_cell.at[0, 'sh_t2id']*1e-3
        sh_t3id = single_cell.at[0, 'sh_t3id']*1e-3
        sh_t4id = single_cell.at[0, 'sh_t4id']*1e-3

        # for 1xop area doesnt matter, tau fit is always good, but exp res from hjcfit is more accurate
        op_p1id = single_cell.at[0, 'op_p1id']/100
        op_t1id = single_cell.at[0, 'op_t1id']*1e-3

        sc_tres = single_cell.at[0, 'tres']*1e-6
        sc_scns = list(single_cell.loc[:, 'file_scn'])
        sc_scns = [name if name.endswith('.SCN') else name + '.SCN' for name in sc_scns]

        # print(sc_scns)
        # print(sc_tres)

        ini_sh_taus = [sh_t1id, sh_t2id, sh_t3id, sh_t4id]
        ini_sh_areas = [sh_p1id, sh_p2id, sh_p3id, sh_p4id]

        ini_op_taus = [op_t1id]
        ini_op_areas = [op_p1id]

        # TODO: some printing and comparison to initials

        fit_sh_taus, fit_sh_areas = fit_ekdist_shuts(sc_scns, sc_tres, ini_sh_taus, ini_sh_areas, '#ff0000', cell + '_shuts.png', 1e0)
        fit_op_taus, fit_op_areas = fit_ekdist_openings(sc_scns, sc_tres, ini_op_taus, ini_op_areas , '#0070c0', cell + '_openings.png', 1e-1)

        cell_result = {'cell_id': cell, 'sh_area_ek_1': fit_sh_areas[0], 'sh_area_ek_2': fit_sh_areas[1], 'sh_area_ek_3': fit_sh_areas[2], 'sh_area_ek_4': fit_sh_areas[3],
                       'sh_tau_ek_1': fit_sh_taus[0], 'sh_tau_ek_2': fit_sh_taus[1], 'sh_tau_ek_3': fit_sh_taus[2], 'sh_tau_ek_4': fit_sh_taus[3],
                       'op_area_ek1': fit_op_areas[0], 'op_tau_ek1': fit_op_taus[0]}
        results.append(cell_result)

    results = pd.DataFrame(results)
    results.to_csv("ekdist_results_e153_2025_CFODD_DISTS.csv")

def fit_ekdist_openings(scns, tres, taus, areas, color, plot_name, lim):

    rec = ekrecord.SingleChannelRecord()
    rec.load_SCN_file(scns)
    rec.tres = tres
    print(rec)
    print('')

    #ekplot.plot_xlog_interval_histogram(rec.opint, rec.tres, shut=False)
    #plt.show()

    expPDF = dceq.MultiExponentialPDF(np.asarray(rec.opint), taus=np.asarray(taus), areas=np.asarray(areas))
    theta = expPDF.theta
    #print('Start LogLikelihood =', expPDF.loglik(theta))
    res = minimize(expPDF.loglik, theta, method='Nelder-Mead')
    #print(res)
    expPDF.theta = res.x

    ekplot.plot_xlog_interval_histogram_fit(rec.opint, rec.tres, expPDF.to_plot, res.x, 2, 2, shut=False)

    ax = plt.gca()
    for line in ax.get_lines():
        line.set_color(color)

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(9)
    ax.set_xlabel('')
    ax.set_ylabel('')

    for axis in ['top', 'right']:
        ax.spines[axis].set_linewidth(0)
    ax.tick_params(width=2, which='both')
    ax.set_xlim([3e-5, lim])
    ax.set_ylim([0, 15])

    plt.tight_layout()
    plt.savefig(plot_name, dpi=300)
    #plt.show()

    return expPDF.taus, expPDF.areas


def fit_ekdist_shuts(scns, tres, taus, areas, color, plot_name, lim):

    rec = ekrecord.SingleChannelRecord()
    rec.load_SCN_file(scns)
    rec.tres = tres
    print(rec)
    print('')

    # ekplot.plot_xlog_interval_histogram(rec.shint, rec.tres, shut=True)
    #plt.show()

    expPDF = dceq.MultiExponentialPDF(np.asarray(rec.shint), taus=np.asarray(taus), areas=np.asarray(areas))
    theta = expPDF.theta
    #print('Start LogLikelihood =', expPDF.loglik(theta))
    res = minimize(expPDF.loglik, theta, method='Nelder-Mead')
    #print(res)
    expPDF.theta = res.x

    ekplot.plot_xlog_interval_histogram_fit(rec.shint, rec.tres, expPDF.to_plot, res.x, 2, 2, shut=True)

    ax = plt.gca()
    for line in ax.get_lines():
        line.set_color(color)

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(9)
    ax.set_xlabel('')
    ax.set_ylabel('')

    for axis in ['top', 'right']:
        ax.spines[axis].set_linewidth(0)
    ax.tick_params(width=2, which='both')
    ax.set_xlim([3e-5, lim])
    ax.set_ylim([0, 25])

    plt.tight_layout()
    plt.savefig(plot_name, dpi=300)
    #plt.show()

    return expPDF.taus, expPDF.areas

fit_from_config_4sh1op("ekdist_config_e153_2025_CFODD_DISTS.csv")

#fit_ekdist_shuts(['53WT3_K1.SCN'], 40e-6, [0.02e-3, 0.156e-3, 0.55e-3, 6.04e-3], [0.67, 0.25, 0.07, 0.006], '#ff0000', 'V53WT3_shuts.png', 1e-1)
#fit_ekdist_openings(['53WT3_K1.SCN'], 40e-6, [1.28e-3, 2.19e-3], [0.48, 0.52], '#0070c0', 'V53WT3_openings.png', 1e-1)

#fit_ekdist_shuts(['kom1_1k.SCN', 'KOM1_2k.SCN', 'KOM1_3k.SCN', 'KOM1_4k.SCN', 'KOM1_5k.SCN', 'KOM1_6k.SCN'], 50e-6, [0.06e-3, 0.65e-3, 3.72e-3, 28.58e-3], [0.58, 0.26, 0.13, 0.02], '#ff0000', 'E270K1_shuts.png', 1e0)
#fit_ekdist_openings(['kom1_1k.SCN', 'KOM1_2k.SCN', 'KOM1_3k.SCN', 'KOM1_4k.SCN', 'KOM1_5k.SCN', 'KOM1_6k.SCN'], 50e-6, [0.82e-3, 4.99e-3], [0.39, 0.61], '#0070c0', 'E270K1_openings.png', 1e-1)

#fit_ekdist_shuts(['KOM2(1).SCN', 'KOM2(2).SCN', 'KOM2(3).SCN', 'KOM2(4).SCN', 'KOM2(5).SCN'], 40e-6, [0.08e-3, 0.37e-3, 2.97e-3, 20.18e-3], [0.53, 0.39, 0.06, 0.02], '#ff0000', 'E270Q1_shuts.png', 1e0)
#fit_ekdist_openings(['KOM2(1).SCN', 'KOM2(2).SCN', 'KOM2(3).SCN', 'KOM2(4).SCN', 'KOM2(5).SCN'], 40e-6, [0.31e-3, 1.67e-3], [0.13, 0.87], '#0070c0', 'E270Q1_openings.png', 1e-1)
#fit_ekdist_openings(['KOM2(1).SCN', 'KOM2(2).SCN', 'KOM2(3).SCN', 'KOM2(5).SCN'], 40e-6, [0.31e-3, 1.67e-3], [0.13, 0.87], '#0070c0', 'E270Q1_openings.png', 1e-1)

#fit_ekdist_shuts(['270F1_K1.SCN', '270F1_K7.SCN'], 35e-6, [0.05e-3, 0.40e-3, 1.73e-3, 16.09e-3], [0.25, 0.57, 0.16, 0.02], '#ff0000', 'E270F1_shuts.png', 1e0)
#fit_ekdist_openings(['270F1_K1.SCN', '270F1_K7.SCN'], 35e-6, [1.38e-3, 5.09e-3], [0.9, 0.1], '#0070c0', 'E270F1_openings.png', 1e-1)

#fit_ekdist_shuts(['2406_C1.SCN', '2406_C2.SCN', '2406_C3.SCN'], 70e-6, [0.25e-3, 2.08e-3, 9.47e-3], [0.26, 0.68, 0.06], '#ff0000', 'E270C1_shuts.png', 1e0)
#fit_ekdist_openings(['2406_C1.SCN', '2406_C2.SCN', '2406_C3.SCN'], 70e-6, [0.16e-3, 0.66e-3], [0.78, 0.22], '#0070c0', 'E270C1_openings.png', 1e-1)

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