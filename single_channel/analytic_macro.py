
# ?

import pandas as pd
import plotly.express as px
from scalcs import mechanism
from scalcs import cjumps
import argparse


class Model:

    def __init__(self, topology, rates):
        """

        :param rates:
        """
        self.model_mechanism = (topology, rates)

    @property
    def model_mechanism(self):
        return self._model_mechanism

    @model_mechanism.setter
    def model_mechanism(self, params):

        if params[0] == 'RAAFOODD':
            self._model_mechanism = self.mechanism_RAAFOODD(params[1])

        elif params[0] == 'RAAFODD':
            self._model_mechanism = self.mechanism_RAAFODD(params[1])

        elif params[0] == 'RAAFOD':
            self._model_mechanism = self.mechanism_RAAFOD(params[1])

    def mechanism_RAAFOD(self, rates):

        mectitle = 'R-A-A2-F-0=D'
        ratetitle = 'starters'

        A2D = mechanism.State('B', 'A2D', 0.0)
        A2O = mechanism.State('A', 'A2O', 50e-12)
        A2F = mechanism.State('B', 'A2F', 0.0)
        A2R = mechanism.State('B', 'A2R', 0.0)
        AR = mechanism.State('B', 'AR', 0.0)
        R = mechanism.State('C', 'R', 0.0)

        RateList = [

            mechanism.Rate(rates['res'], A2D, A2F, name='res', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates['des'], A2F, A2D, name='des', limits=[1e-15, 1e+7]),

            mechanism.Rate(rates['alpha'], A2O, A2F, name='alpha', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates['beta'], A2F, A2O, name='beta', limits=[1e-15, 1e+7]),

            mechanism.Rate(rates['gamma'], A2F, A2R, name='gamma', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates['delta'], A2R, A2F, name='delta', limits=[1e-15, 1e+7]),

            mechanism.Rate(rates['koff'] * 2, A2R, AR, name='2koff', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates['kon'] * 2, R, AR, name='2kon', eff='c', limits=[1e-15, 1e+10]),
            mechanism.Rate(rates['koff'], AR, R, name='koff', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates['kon'], AR, A2R, name='kon', eff='c', limits=[1e-15, 1e+10]),
        ]

        return mechanism.Mechanism(RateList, mtitle=mectitle, rtitle=ratetitle)

    def mechanism_RAAFODD(self, rates):

        mectitle = 'R-A-A2-F-0=2D'
        ratetitle = 'starters'

        A2Dp = mechanism.State('B', 'A2Dp', 0.0)
        A2D = mechanism.State('B', 'A2D', 0.0)
        A2O = mechanism.State('A', 'A2O', 50e-12)
        A2F = mechanism.State('B', 'A2F', 0.0)
        A2R = mechanism.State('B', 'A2R', 0.0)
        AR = mechanism.State('B', 'AR', 0.0)
        R = mechanism.State('C', 'R', 0.0)

        RateList = [

            mechanism.Rate(rates['res'], A2D, A2F, name='res', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates['des'], A2F, A2D, name='des', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates['res_p'], A2Dp, A2F, name='res_p', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates['des_p'], A2F, A2Dp, name='des_p', limits=[1e-15, 1e+7]),

            mechanism.Rate(rates['alpha'], A2O, A2F, name='alpha', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates['beta'], A2F, A2O, name='beta', limits=[1e-15, 1e+7]),

            mechanism.Rate(rates['gamma'], A2F, A2R, name='gamma', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates['delta'], A2R, A2F, name='delta', limits=[1e-15, 1e+7]),

            mechanism.Rate(rates['koff'] * 2, A2R, AR, name='2koff', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates['kon'] * 2, R, AR, name='2kon', eff='c', limits=[1e-15, 1e+10]),
            mechanism.Rate(rates['koff'], AR, R, name='koff', limits=[1e-15, 1e+7]),
            mechanism.Rate(rates['kon'], AR, A2R, name='kon', eff='c', limits=[1e-15, 1e+10]),
        ]

        return mechanism.Mechanism(RateList, mtitle=mectitle, rtitle=ratetitle)

    def mechanism_RAAFOODD(self, rates):

        mectitle = 'R-A-A2-F-20=2D'
        ratetitle = 'starters'

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

        return mechanism.Mechanism(RateList, mtitle=mectitle, rtitle=ratetitle)


class ModelsBuilder:

    # def __init__(self, s_rates, v_rate, topo):
    def __init__(self, s_rates, topo):

        self.topology = topo
        self.start_rates = s_rates
        # self.variable_rate = v_rate

    @property
    def topology(self):
        return self._topology

    @topology.setter
    def topology(self, topo):
        self._topology = topo

    @property
    def start_rates(self):
        return self._start_rates

    @start_rates.setter
    def start_rates(self, s_rates):
        self._start_rates = s_rates

    '''
    @property
    def variable_rate(self):
        return self._variable_rate

    @variable_rate.setter
    def variable_rate(self, v_rate):
        self._variable_rate = v_rate
    '''

    '''
    def build_models(self):

        step_size = 8e-5
        record_length = 1000e-3

        model_trace = pd.DataFrame()

        for variable_rate_val in [0.1 * self.start_rates[self.variable_rate], 0.3 * self.start_rates[self.variable_rate],
                                  0.5 * self.start_rates[self.variable_rate], 0.7 * self.start_rates[self.variable_rate],
                                  0.9 * self.start_rates[self.variable_rate],
                                  self.start_rates[self.variable_rate],
                                  1.1 * self.start_rates[self.variable_rate], 1.3 * self.start_rates[self.variable_rate],
                                  1.5 * self.start_rates[self.variable_rate], 1.7 * self.start_rates[self.variable_rate],
                                  2.0 * self.start_rates[self.variable_rate]]:

            self.start_rates[self.variable_rate] = variable_rate_val
            sample_model = Model(self.topology, start_rates)

            t, c, Popen, P = cjumps.solve_jump(sample_model.model_mechanism, record_length, step_size, cjumps.pulse_square,
                                               (1e-2, 0.0, 100e-3, 500e-3))

            model_trace[variable_rate_val] = Popen

        model_trace['t'] = t * 1000
        model_trace.set_index(keys='t', drop=True, inplace=True)

        model_trace_normalized = model_trace.divide(model_trace.max(), axis=1)
        model_trace_normalized_forATF = model_trace_normalized.copy()
        print(model_trace_normalized_forATF)

        model_trace_normalized.reset_index(inplace=True)
        model_trace_normalized_long = (
            model_trace_normalized.melt(id_vars=['t'], var_name=variable_rate, value_name='Popen'))

        fig = px.line(model_trace_normalized_long, x='t', y='Popen', color=variable_rate,
                      color_discrete_sequence=px.colors.qualitative.Dark24, title=str(start_rates))
        fig.write_html(self.variable_rate + '.html')

        self.save_atf('header.txt', self.variable_rate + '.atf', model_trace_normalized_forATF)
    '''

    def build_models_multi(self):

        step_size = 8e-5
        record_length = 1000e-3
        model_traces = []

        for variable_rate in self.start_rates.keys():

            steps = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 2.0, 5.0]

            single_var_trace_atf = pd.DataFrame()

            for step in steps:

                variable_rate_val = self.start_rates[variable_rate] * step

                current_rates = self.start_rates.copy()
                current_rates[variable_rate] = variable_rate_val
                sample_model = Model(self.topology, current_rates)

                t, c, Popen, P = cjumps.solve_jump(sample_model.model_mechanism, record_length, step_size,
                                                   cjumps.pulse_square, (1e-2, 0.0, 100e-3, 500e-3))

                model_trace = pd.DataFrame()

                model_trace['Popen'] = Popen
                model_trace = model_trace.divide(model_trace.max(), axis=1)
                model_trace['t'] = t * 1000

                model_trace['Rate_name'] = variable_rate
                model_trace['Rate_value'] = variable_rate_val
                model_trace['Rate_value_step'] = step

                print(variable_rate, variable_rate_val)
                print(model_trace)
                model_traces.append(model_trace)

                single_var_trace_atf[variable_rate_val] = Popen

            single_var_trace_atf['t'] = t * 1000
            single_var_trace_atf.set_index(keys='t', drop=True, inplace=True)
            self.save_atf('header.txt', variable_rate + '.atf', single_var_trace_atf)

        all_traces = pd.concat(model_traces)

        fig = px.line(all_traces, x='t', y='Popen', color='Rate_value_step', facet_col='Rate_name', facet_col_wrap=2,
                      height=500 * len(self.start_rates.keys()),
                      template='simple_white',
                      color_discrete_sequence=px.colors.sequential.RdBu,
                      hover_data=['Rate_value'])

        fig.update_layout(font=dict(family="Courier New, monospace", size=18))

        fig.write_html('simulations' + '.html')

    def save_atf(self, header, outfile, multi_data):
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


parser = argparse.ArgumentParser()
# parser.add_argument("-vr", "--variable_rate", type=str)
parser.add_argument("-sr", "--start_rates", type=str)
parser.add_argument("-tp", "--topology", type=str)


args = parser.parse_args()


start_rates = pd.read_csv(args.start_rates, header=None)
start_rates = dict(zip(start_rates[0], start_rates[1]))
# variable_rate = args.variable_rate

# builder = ModelsBuilder(start_rates, variable_rate, args.topology)

builder = ModelsBuilder(start_rates, args.topology)
#builder.build_models()
builder.build_models_multi()




