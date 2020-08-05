
# ?

import pandas as pd
import plotly.express as px
from scalcs import mechanism
from scalcs import cjumps


class ModelCCOD:

    def __init__(self, rates):
        """

        :param rates:
        """

        self.mechanism = rates

    @property
    def mechanism(self):
        return self._mechanism

    @mechanism.setter
    def mechanism(self, rates):

        mectitle = 'C-C-O-D'
        ratetitle = 'quasi random numbers'

        A2Dp = mechanism.State('B', 'A2Dp', 0.0)
        A2D  = mechanism.State('B', 'A2D', 0.0)
        A2Op = mechanism.State('A', 'A2Op', 50e-12)
        A2O  = mechanism.State('A', 'A2O', 50e-12)
        A2F  = mechanism.State('B', 'A2F', 0.0)
        A2R  = mechanism.State('B', 'A2R', 0.0)
        AR   = mechanism.State('B', 'AR', 0.0)
        R    = mechanism.State('C', 'R', 0.0)

        RateList = [

             mechanism.Rate(rates['res'], A2D, A2F, name='res', limits=[1e-15, 1e+7]),
             mechanism.Rate(rates['des'], A2F, A2D, name='des', limits=[1e-15, 1e+7]),
             mechanism.Rate(rates['res_p'], A2Dp, A2F, name='res_p', limits=[1e-15, 1e+7]),
             mechanism.Rate(rates['des_p'], A2F, A2Dp, name='des_p', limits=[1e-15, 1e+7]),

             mechanism.Rate(rates['alpha'], A2O, A2F, name='alpha', limits=[1e-15,1e+7]),
             mechanism.Rate(rates['beta'], A2F, A2O, name='beta', limits=[1e-15,1e+7]),
             mechanism.Rate(rates['alpha_p'], A2Op, A2F, name='alpha_p', limits=[1e-15, 1e+7]),
             mechanism.Rate(rates['beta_p'], A2F, A2Op, name='beta_p', limits=[1e-15, 1e+7]),

             mechanism.Rate(rates['gamma'], A2F, A2R, name='gamma', limits=[1e-15,1e+7]),
             mechanism.Rate(rates['delta'], A2R, A2F, name='delta', limits=[1e-15,1e+7]),

             mechanism.Rate(rates['koff']*2, A2R, AR, name='2koff', limits=[1e-15, 1e+7]),
             mechanism.Rate(rates['kon']*2, R, AR, name='2kon', eff='c', limits=[1e-15, 1e+10]),
             mechanism.Rate(rates['koff'], AR, R, name='koff', limits=[1e-15,1e+7]),
             mechanism.Rate(rates['kon'], AR, A2R, name='kon', eff='c', limits=[1e-15,1e+10]),
             ]

        self._mechanism = mechanism.Mechanism(RateList, mtitle=mectitle, rtitle=ratetitle)


def to_pc(header, outfile, multi_data):
    '''
    saves Pandas data frame as single atf file
    :param header: path to header template
    :param outfile: path to atf outfile
    :param multi_data: Pandas data frame with all the abfs data
    :return:
    '''

    # prepare header
    with open(header, 'r') as h:
        cols = multi_data.shape[1]

        header = h.readlines()
        header[1] = "7\t{}\n".format(1 + cols)
        header[8] = header[8].rstrip() + "\t" + "\"sim #1\"\t" * cols + "\n"
        header[9] = header[9].rstrip() + "\t" + "\"Trace #1 (prob)\"\t" * cols + "\n"
        header = ''.join(header)

    # write header
    with open(outfile, 'w') as d:
        d.write(header)

    # write all data
    with open(outfile, 'a') as d:
        multi_data.to_csv(d, header=False, float_format='%.4f', sep='\t')


start_rates = {'kon': 10e6, 'koff': 900,
               'delta': 4400, 'gamma': 3000,
               'beta': 6900, 'alpha': 1800, 'beta_p': 4300, 'alpha_p': 450,
               'des': 4000, 'res': 60, 'des_p': 400, 'res_p': 3}

step_size = 8e-5
record_length = 1000e-3
model_trace = pd.DataFrame()

variable_rate = 'delta'

for variable_rate_val in [0.1*start_rates[variable_rate], 0.3*start_rates[variable_rate], 0.5*start_rates[variable_rate], 0.7*start_rates[variable_rate], 0.9*start_rates[variable_rate],
                          start_rates[variable_rate],
                          1.1*start_rates[variable_rate], 1.3*start_rates[variable_rate], 1.5*start_rates[variable_rate], 1.7*start_rates[variable_rate], 2.0*start_rates[variable_rate]]:

    start_rates[variable_rate] = variable_rate_val
    sample_model = ModelCCOD(start_rates)

    t, c, Popen, P = cjumps.solve_jump(sample_model.mechanism, record_length, step_size, cjumps.pulse_square, (1e-2, 0.0, 100e-3, 500e-3))

    model_trace[variable_rate_val] = Popen

model_trace['t'] = t * 1000
model_trace.set_index(keys='t', drop=True, inplace=True)

model_trace_normalized = model_trace.divide(model_trace.max(), axis=1)
model_trace_normalized_forATF = model_trace_normalized.copy()


model_trace_normalized.reset_index(inplace=True)
model_trace_normalized_long = (model_trace_normalized.melt(id_vars=['t'], var_name=variable_rate, value_name='Popen'))


print(model_trace_normalized_forATF)

fig = px.line(model_trace_normalized_long, x='t', y='Popen', color=variable_rate, color_discrete_sequence=px.colors.qualitative.Dark24, title=str(start_rates))
fig.write_html('test.html')

to_pc('header.txt', 'test.atf', model_trace_normalized_forATF)