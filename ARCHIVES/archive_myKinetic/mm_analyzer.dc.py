import logging as log
import numpy as np

from dcpyps import dcio
from dcpyps import dataset
from dcpyps import mechanism
from dcpyps import fitMLL


class KineticDc:

    def __init__(self, mec_file_name, mec_model_num, fixed, scn_files, scn_tres, scn_tcrit, scn_conc):

        # mechanism preparation

        self.mec_file_name, self.mec_model_num = mec_file_name, mec_model_num
        self.version, self.mec_list, self.mec_num = dcio.mec_get_list(self.mec_file_name)

        self.mec = dcio.mec_load(self.mec_file_name, self.mec_list[0][self.mec_model_num])
        for rate in self.mec.Rates:
            if rate._get_name().strip() in fixed:
                log.info('Fixing rate: {}'.format(rate._get_name()))
                rate.fixed = True
            else:
                log.info('Not fixing rate: {}'.format(rate._get_name()))
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
        log.info('Mec info: {}'.format(self.mec))

        # data preparation

        self.scn_files, self.scn_tres, self.scn_tcrit, self.scn_conc = scn_files, scn_tres, scn_tcrit, scn_conc
        self.recs = []
        self.bursts = []

        for sfs, cn, tr, tc in zip(self.scn_files, self.scn_conc, self.scn_tres, self.scn_tcrit):
            rec = dataset.SCRecord(sfs, cn, tr, tc)
            rec.record_type = 'recorded'
            self.recs.append(rec)
            self.bursts.append(rec.bursts.intervals())
            log.info(rec)

        log.info("Reading files: {} with t_res: {} s, t_crit {} s and concentrations: {} M".format(self.scn_files, self.scn_tres, self.scn_tcrit, self.scn_conc))

log.basicConfig(filename='mm_dc.log', filemode='w', format='%(message)s', level=log.DEBUG)
model_dc = KineticDc('magda.mec', 0,
                  [],
                  [['121015C8.SCN', '12101591.SCN', '12101592.SCN', '13101511.SCN', '13101521.SCN', '13101531.SCN'],
                   ['103161K1.SCN', '103161K2.SCN', '103161S1.SCN', '103161S2.SCN', '103162K1.SCN'], #, '103162S1.SCN', '103162S2.SCN'],
                   ['16715151.SCN', '16715152.SCN'],
                   ['18915151.SCN', '18915152.SCN', '18915153.SCN', '18915154.SCN', '18915155.SCN', '18915156.SCN']],
                  [0.000065, 0.000065, 0.000065, 0.000065],
                  [0.0105, 0.007, 0.005, 0.006],
                  [10e-3, 10e-6, 1e-6, 0.1e-6])

fit_dc = fitMLL.FittingSession(model_dc.mec, model_dc.recs)
fit_dc.run()
