import logging as log
import numpy as np
import time
import math
from scipy.optimize import minimize


from dcpyps import dcio
from dcpyps import dataset
from dcpyps import mechanism
from dcprogs.likelihood import Log10Likelihood


class fitterDC:

    def __init__(self, mec_file_name, mec_model_num, fixed, scn_files, scn_tres, scn_tcrit, scn_conc):

        # mechanism preparation

        self.mec_file_name, self.mec_model_num = mec_file_name, mec_model_num
        self.version, self.mec_list, self.mec_num = dcio.mec_get_list(self.mec_file_name)

        self.mec = dcio.mec_load(self.mec_file_name, self.mec_list[0][self.mec_model_num])
        for rate in self.mec.Rates:
            if rate._get_name().strip() in fixed:
                print('Fixing rate: {}'.format(rate._get_name()))
                rate.fixed = True
            else:
                print('Not fixing rate: {}'.format(rate._get_name()))
                rate.fixed = False

        self.mec.Rates[14].is_constrained = True
        self.mec.Rates[14].constrain_func = mechanism.constrain_rate_multiple
        self.mec.Rates[14].constrain_args = [16, 2]

        self.mec.Rates[17].is_constrained = True
        self.mec.Rates[17].constrain_func = mechanism.constrain_rate_multiple
        self.mec.Rates[17].constrain_args = [15, 2]

        self.theta = np.log(self.mec.theta())

        log.info('Mec file version: {}'.format(self.version))
        log.info('Mec list (start byte, model number and description): {}'.format(self.mec_list))
        log.info('Total number of models: {}'.format(self.mec_num))
        log.info('Mec info: {}'.format(self.mec.printout()))

        # data preparation

        self.scn_files, self.scn_tres, self.scn_tcrit, self.scn_conc = scn_files, scn_tres, scn_tcrit, scn_conc
        self.recs = []
        self.bursts = []

        for sfs, cn, tr, tc in zip(self.scn_files, self.scn_conc, self.scn_tres, self.scn_tcrit):
            rec = dataset.SCRecord(sfs, cn, tr, tc)
            rec.record_type = 'recorded'
            self.recs.append(rec)
            self.bursts.append(rec.bursts.intervals())
            log.info(rec.printout())

        log.info("Reading files: {} with t_res: {} s, t_crit {} s and concentrations: {} M".format(self.scn_files, self.scn_tres, self.scn_tcrit, self.scn_conc))

    def optimize_rates(self):

        kwargs = {'nmax': 2, 'xtol': 1e-12, 'rtol': 1e-12, 'itermax': 100, 'lower_bound': -1e6, 'upper_bound': 0}

        likelihood = []

        for rec, bur in zip(self.recs, self.bursts):
            likelihood.append(Log10Likelihood(bur, self.mec.kA, rec.tres, rec.tcrit, **kwargs))

        def dcprogslik(x, args=None):
            self.mec.theta_unsqueeze(np.exp(x))
            lik = 0
            for idx, con in enumerate(self.scn_conc):
                self.mec.set_eff('c', con)
                lik += -likelihood[idx](self.mec.Q) * math.log(10)
            return lik

        self.iternum = 0

        def printiter(theta):
            self.iternum += 1
            lik = dcprogslik(theta)
            print("iteration # {0:d}; log-lik = {1:.6f}".format(self.iternum, -lik))
            print(np.exp(theta))

        lik = dcprogslik(self.theta)
        print("\nStarting likelihood (DCprogs)= {0:.6f}".format(-lik))
        start = time.clock()
        success = False

        while not success:
            result = minimize(dcprogslik, self.theta, method='Nelder-Mead', callback=printiter, options={'xtol':1e-4, 'ftol':1e-4, 'maxiter': 5000, 'maxfev': 10000, 'disp': True})
        if result.success:
            success = True
        else:
            self.theta = result.x

        end = time.clock()
        print('time in simplex='.format(end - start))
        print('\n\nresult=')
        print(result)

        print('\n Final log-likelihood = {0:.6f}'.format(-result.fun))
        print('\n Number of iterations = {0:d}'.format(result.nit))
        print('\n Number of evaluations = {0:d}'.format(result.nfev))
        self.mec.theta_unsqueeze(np.exp(result.x))
        print("\n Final rate constants:")
        self.mec.printout()
        print('\n\n')



log.basicConfig(format='%(message)s', level=log.DEBUG)
#fitter = fitterDC('magda.mec', 0,
#                  ['alfa1', 'alfa2', 'alfa3', 'alfa4', 'beta1', 'beta2', 'beta3', 'beta4', 'k on', 'k off', 'd1', 'r1', 'd2', 'r2'],
#                  [['121015C8.SCN'], ['18915151.SCN']],
#                  [0.000065, 0.000065], [0.0105, 0.006], [10e-3, 0.1e-6])

fitter = fitterDC('magda.mec', 0,
                  ['alfa1', 'alfa2', 'alfa3', 'alfa4', 'beta1', 'beta2', 'beta3', 'beta4', 'k on', 'k off', 'd1', 'r1', 'd2', 'r2'],
                  [['121015C8.SCN']],
                  [0.000065], [0.0105], [10e-3])

fitter.optimize_rates()
