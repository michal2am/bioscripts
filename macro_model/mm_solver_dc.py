import logging as log
import numpy as np
import time
import math
from scipy.optimize import minimize


from dcpyps import dcio
from dcpyps import dataset
from dcpyps import mechanism
from dcprogs.likelihood import Log10Likelihood

log.basicConfig(format='%(message)s', level=log.DEBUG)

# model read

mec_file_name = 'magda.mec'
version, mec_list, mec_num = dcio.mec_get_list(mec_file_name)

log.info('Mec file version: {}'.format(version))
log.info('Mec list (start byte, model number and description): {}'.format(mec_list))
log.info('Total number of models: {}'.format(mec_num))

selected_model = 0
start_byte_index = 0

mec = dcio.mec_load(mec_file_name, mec_list[start_byte_index][selected_model])

for rate in mec.Rates:
    rate.fixed = False

theta = np.log(mec.theta())

log.info('Mec info: {}'.format(mec.printout()))

# trajectory read

scnfiles = [['121015C8.SCN'], ['18915151.SCN']]
tres = [0.000065, 0.000065]
tcrit = [0.0105, 0.006]
conc = [10e-3, 0.1e-6]

log.info("Reading files: {} with t_res: {} s, t_crit {} s and concentrations: {} M".format(scnfiles, tres, tcrit, conc))

recs = []
bursts = []

for sfs, cn, tr, tc in zip(scnfiles, conc, tres, tcrit):
    rec = dataset.SCRecord(sfs, cn, tr, tc)
    rec.record_type = 'recorded'
    recs.append(rec)
    bursts.append(rec.bursts.intervals())
    log.info(rec.printout())

# solving

kwargs = {'nmax': 2, 'xtol': 1e-12, 'rtol': 1e-12, 'itermax': 100, 'lower_bound': -1e6, 'upper_bound': 0}

likelihood = []

for i in range(len(recs)):
    likelihood.append(Log10Likelihood(bursts[i], mec.kA, recs[i].tres, recs[i].tcrit, **kwargs))

def dcprogslik(x, args=None):
    mec.theta_unsqueeze(np.exp(x))
    lik = 0
    for i in range(len(conc)):
        mec.set_eff('c', conc[i])
        lik += -likelihood[i](mec.Q) * math.log(10)
    return lik

iternum = 0
def printiter(theta):
    global iternum
    iternum += 1
    lik = dcprogslik(theta)
    print("iteration # {0:d}; log-lik = {1:.6f}".format(iternum, -lik))
    print(np.exp(theta))

lik = dcprogslik(theta)
print("\nStarting likelihood (DCprogs)= {0:.6f}".format(-lik))
start = time.clock()
success = False

while not success:
    result = minimize(dcprogslik, theta, method='Nelder-Mead', callback=printiter, options={'xtol':1e-4, 'ftol':1e-4, 'maxiter': 5000, 'maxfev': 10000, 'disp': True})
