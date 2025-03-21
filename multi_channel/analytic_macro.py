# performs automated macro simulations + parameter estimation
# example call: analytic_macro.py --start_rates start_rates_RAAFOODD_makIlo_schrels_d396dp185.csv --topology RAAFOODD --annotation wt_makIlo_schrels_d396dp185

import pandas as pd
import numpy as np
import plotly.express as px
# import matplotlib.pyplot as plt
import argparse
from scalcs import mechanism
from scalcs import cjumps
from scalcs import popen
from scipy.optimize import curve_fit


class Model:

    def __init__(self, topology, rates):
        """
        builds a single SCALCS model
        str:param topology: RAAFOD, RAAFODD or RAAFOODD:
        dict:param rates: compatible with selected topology
        """
        self.model_mechanism = (topology, rates)

    @property
    def model_mechanism(self):
        """
        :return: SCALCS model
        """
        return self._model_mechanism

    @model_mechanism.setter
    def model_mechanism(self, params):
        """"
        list:param params: topology and rates
        """

        if params[0] == 'RAAFOODD':
            self._model_mechanism = self.mechanism_RAAFOODD(params[1])

        elif params[0] == 'RAAFODD':
            self._model_mechanism = self.mechanism_RAAFODD(params[1])

        elif params[0] == 'RAAFOD':
            self._model_mechanism = self.mechanism_RAAFOD(params[1])

    @staticmethod
    def mechanism_RAAFOD(rates):
        """
        dict:param rates: rates
        :return: minimal SCALCS model
        """

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

    @staticmethod
    def mechanism_RAAFODD(rates):
        """
        dict:param rates: rates
        :return: double desensitization SCALCS model
        """

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

    @staticmethod
    def mechanism_RAAFOODD(rates):
        """
        dict:param rates: rates
        :return: double desensitization and open SCALCS model
        """

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

             mechanism.Rate(rates['alpha'], A2O, A2F, name='alpha', limits=[1e-15, 1e+7]),
             mechanism.Rate(rates['beta'], A2F, A2O, name='beta', limits=[1e-15, 1e+7]),
             mechanism.Rate(rates['alpha_p'], A2Op, A2F, name='alpha_p', limits=[1e-15, 1e+7]),
             mechanism.Rate(rates['beta_p'], A2F, A2Op, name='beta_p', limits=[1e-15, 1e+7]),

             mechanism.Rate(rates['gamma'], A2F, A2R, name='gamma', limits=[1e-15, 1e+7]),
             mechanism.Rate(rates['delta'], A2R, A2F, name='delta', limits=[1e-15, 1e+7]),

             mechanism.Rate(rates['koff']*2, A2R, AR, name='2koff', limits=[1e-15, 1e+7]),
             mechanism.Rate(rates['kon']*2, R, AR, name='2kon', eff='c', limits=[1e-15, 1e+10]),
             mechanism.Rate(rates['koff'], AR, R, name='koff', limits=[1e-15, 1e+7]),
             mechanism.Rate(rates['kon'], AR, A2R, name='kon', eff='c', limits=[1e-15, 1e+10]),
             ]

        return mechanism.Mechanism(RateList, mtitle=mectitle, rtitle=ratetitle)


class ModelsBuilder:

    def __init__(self, s_rates, e_rates, topo, annotation, concentration, pulse_length):
        """
        builds multiple SCALCS models and performs basic analysis/simulations
        dict:param s_rates: starting rates
        list:param e_rates: rates to exlude from mutli value simulations
        str:param topo: model topology
        """
        self.topology = topo
        self.start_rates = s_rates
        self.exclude_rates = e_rates
        self.annotation = annotation
        self.concentration = concentration
        self.pulse_length = pulse_length

    # TODO: add missing get/set property methods

    @property
    def topology(self):
        """
        str:return: topology
        """
        return self._topology

    @topology.setter
    def topology(self, topo):
        """
        str:param topo: topology
        """
        self._topology = topo

    @property
    def start_rates(self):
        """
        dict:return: starting rates
        """
        return self._start_rates

    @start_rates.setter
    def start_rates(self, s_rates):
        """
        dict:param s_rates: starting rates
        """
        self._start_rates = s_rates

    @staticmethod
    def fit_rise(trace_rise):
        """
        dataframe:param trace_rise: cut and 0-started rise period
        float:return: rise time
        """

        def func_rise(x, a):
            return 1 - np.exp(-x / a)

        try:
            rise_popt, rise_pcov = curve_fit(func_rise, trace_rise['t'], trace_rise['Popen'], p0=[0.25])
        except RuntimeError:
            print('FIT ERROR => RT = 0')
            rise_popt = [0]

        rt = rise_popt[0] * np.log(9)
        print('Rise time:')
        print(rt)

        return '{:.2f}'.format(rt)

    @staticmethod
    def fit_des(trace_des):
        """
        dataframe:param trace_des: cut and 0-started desensitization period
        list:return: desensitization parameters
        """

        def func_des(x, a, b, c, d, e):
            return a * np.exp(-x / b) + c * np.exp(-x / d) + e

        try:
            des_popt, des_pcov = curve_fit(func_des, trace_des['t'], trace_des['Popen'],
                                           p0=(0.75, 2.0, 0.25, 150.0, 0.2))
        except RuntimeError:
            print('FIT ERROR => desensitization parameters = 0')
            des_popt = [0, 0, 0, 0, 0]

        d_a1, d_t1, d_a2, d_t2, d_a3 = des_popt
        print('Desensitization parameters, A1, t1, A2, t2, A3:')
        print(des_popt)

        return ['{:.2f}'.format(param) for param in [d_a1, d_t1, d_a2, d_t2, d_a3]]

    @staticmethod
    def fit_des_single(trace_des):
        """
        dataframe:param trace_des: cut and 0-started desensitization period (single phase)
        list:return: desensitization parameters
        """

        def func_des(x, a, b, c):
            return a * np.exp(-x / b) + c

        try:
            des_popt, des_pcov = curve_fit(func_des, trace_des['t'], trace_des['Popen'], p0=(1, 2.0, 0.2))
        except RuntimeError:
            print('FIT ERROR => desensitization parameters = 0')
            des_popt = [0, 0, 0]

        d_a1, d_t1, d_a3 = des_popt
        print('Desensitization parameters, A1, t1, A2, t2, A3:')
        print(des_popt)

        return ['{:.2f}'.format(param) for param in [d_a1, d_t1, 0, 0, d_a3]]

    @staticmethod
    def fit_dea(trace_dea):
        """
        dataframe:param trace_dea: cut and 0-started deactivation period
        float:return: deactivation time constant
        """

        def func_dea(x, a, b, c, d):
            return a * np.exp(-x / b) + c * np.exp(-x / d)

        try:
            dea_popt, dea_pcov = curve_fit(func_dea, trace_dea['t'], trace_dea['Popen'], p0=[5e-2, 5e2, 6e-2, 7e1])
        except RuntimeError:
            print('FIT ERROR => deactivation parameters = 0')
            dea_popt = [1, 0, 1, 0]

        print('Deactivation parameters, A1, t1, A2, t2:')
        print(dea_popt)
        dea_mean = dea_popt[1] * (dea_popt[0] / (dea_popt[0] + dea_popt[2])) + dea_popt[3] * (
                    dea_popt[2] * (dea_popt[0] + dea_popt[2]))
        print('Deactivation mean:')
        print(dea_mean)

        return '{:.2f}'.format(dea_mean)

    @staticmethod
    def save_atf(header, outfile, multi_data):
        """
        saves Pandas data frame as single atf file
        :param header: path to header template
        :param outfile: path to atf outfile
        :param multi_data: Pandas data frame with all the abfs data
        """

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

    def build_models_multi(self):
        """
        builds multiple SCALCS model for macroscopic trend analysis
        """
        step_size = 8e-5
        # TODO: parametrize record_length
        record_length = self.pulse_length + 2000e-3     # add 1000 ms after pulse
        model_traces = []

        for variable_rate in self.start_rates.keys():

            if self.exclude_rates and variable_rate in self.exclude_rates:
                steps = [1]
            else:
                steps = [0.05, 0.10, 0.15, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4]
                #steps = [0.1, 0.25, 0.5, 0.75, 1, 1.3, 2, 4, 10]
                #steps = [1]

            single_var_trace_atf = pd.DataFrame()

            for step in steps:

                variable_rate_val = self.start_rates[variable_rate] * step

                current_rates = self.start_rates.copy()
                current_rates[variable_rate] = variable_rate_val
                sample_model = Model(self.topology, current_rates)

                # concentration here?
                # problem with lower than 1e-4
                # pulse_square : max concentration / background concentration / prepulse / pulse
                # TODO: parametrize pulse length
                t, c, p_open, p = cjumps.solve_jump(sample_model.model_mechanism, record_length, step_size,
                                                    cjumps.pulse_square, (self.concentration, 0.0, 100e-3,
                                                                          self.pulse_length))

                model_trace = pd.DataFrame()

                model_trace['Popen'] = p_open
                max_a = model_trace.max()
                model_trace = model_trace.divide(max_a, axis=1)
                model_trace['t'] = t * 1000
                model_trace['a'] = max_a.at['Popen']
                
                print('EC50')
                print(popen.EC50(sample_model.model_mechanism, 0))
                model_trace['EC50'] = (popen.EC50(sample_model.model_mechanism, 0))

                model_trace['Rate_name'] = variable_rate
                model_trace['Rate_value'] = variable_rate_val
                model_trace['Rate_value_step'] = step

                print('### Another set ###')
                print('rate and vaule:')
                print(variable_rate, variable_rate_val)
                print('full simulated trace:')
                print(model_trace)

                trace_rise = model_trace.iloc[1250:model_trace['Popen'].idxmax()].copy()
                trace_rise.loc[:, 't'] -= 100.00
                param_rt = self.fit_rise(trace_rise)
                model_trace['rt'] = float(param_rt)

                if self.topology == 'RAAFOD':
                    # conc
                    pass
                    #trace_des = model_trace.iloc[model_trace['Popen'].idxmax():1875].copy()
                    #trace_des.loc[:, 't'] -= trace_des['t'][trace_des['Popen'].idxmax()]
                    #params_des = self.fit_des_single(trace_des)
                else:
                    pass
                    #trace_des = model_trace.iloc[model_trace['Popen'].idxmax():7500].copy()
                    #trace_des.loc[:, 't'] -= trace_des['t'][trace_des['Popen'].idxmax()]
                    #params_des = self.fit_des(trace_des)

                # conc
                #model_trace['d_a1'] = float(params_des[0])
                #model_trace['d_t1'] = float(params_des[1])
                #model_trace['d_a2'] = float(params_des[2])
                #model_trace['d_t2'] = float(params_des[3])
                #model_trace['d_a3'] = float(params_des[4])

                trace_dea = model_trace.iloc[7500:].copy()
                trace_dea.loc[:, 't'] -= 600
                param_dea = self.fit_dea(trace_dea)
                model_trace['dea_m'] = float(param_dea)

                model_traces.append(model_trace)
                single_var_trace_atf[variable_rate_val] = p_open

            single_var_trace_atf['t'] = t * 1000
            single_var_trace_atf.set_index(keys='t', drop=True, inplace=True)
            self.save_atf('header.txt', self.topology + '_' + variable_rate + '.atf', single_var_trace_atf)

        all_traces = pd.concat(model_traces)

        parameters = all_traces.drop_duplicates(subset=['Rate_name', 'Rate_value'])
        parameters = parameters.drop(labels=['Popen', 't'], axis=1)
        parameters = parameters.melt(id_vars=['Rate_name', 'Rate_value', 'Rate_value_step'])

        parameters.to_csv('parameters.csv')

        # Somehow got messed up - y-axis is not ordered, maybe values are not floats?
        # fixed by float() where params are put into the df
        # print(parameters)
        # parameters = pd.to_numeric(parameters['value'], errors='coerce')
        # print(parameters.dtypes)

        fig2 = px.scatter(parameters, x='Rate_value_step', y='value', facet_col='variable', color='Rate_name',
                          facet_col_wrap=4,
                          height=1000,
                          template='simple_white',
                          color_discrete_sequence=px.colors.qualitative.D3,
                          hover_data=['Rate_value']
                          )

        fig2.update_yaxes(matches=None, showticklabels=True)
        fig2.update_traces(mode='lines+markers')
        fig2.update_layout(font=dict(family="Courier New, monospace", size=18))
        fig2.write_html(self.topology + '_' + self.annotation + '_parameters' + '.html')
        # fig2.write_image(self.topology + '_parameters' + '.png')


        fig = px.line(all_traces, x='t', y='Popen', facet_col='Rate_name', color='Rate_value_step',
                      facet_col_wrap=4,
                      height=125 * len(self.start_rates.keys()),
                      template='simple_white',
                      color_discrete_sequence=px.colors.sample_colorscale("amp", [n/(len(steps) -1) for n in range(len(steps))]))
                      # conc
                      # hover_data=['Rate_value', 'rt', 'd_a1', 'd_t1', 'd_a2', 'd_t2', 'd_a3', 'dea_m'],)

        fig.update_layout(font=dict(family="Courier New, monospace", size=18))

        fig.write_html(self.topology + '_' + self.annotation + '_simulations' + '.html')
        # fig.write_image(self.topology + '_simulations' + '.png')


parser = argparse.ArgumentParser()
parser.add_argument("-sr", "--start_rates", type=str)
parser.add_argument("-er", "--exclude_rates", nargs='+')
parser.add_argument("-tp", "--topology", type=str)
parser.add_argument("-an", "--annotation", type=str)
parser.add_argument("-co", "--concentration", type=float)
parser.add_argument("-pl", "--pulse_length", type=float)

args = parser.parse_args()

start_rates = pd.read_csv(args.start_rates, header=None)
start_rates = dict(zip(start_rates[0], start_rates[1]))

builder = ModelsBuilder(start_rates, args.exclude_rates, args.topology, args.annotation, args.concentration,
                        args.pulse_length)
builder.build_models_multi()


'''

plot_fits = False
if plot_fits:
    # plot to check desensitization fit


    # plot to check desensitization fit
    curvey = func_dea(trace_dea['t'], dea_popt[0], dea_popt[1], dea_popt[2],
                      dea_popt[3])  # This is your y axis fit-line
    plt.plot(trace_dea['t'], curvey, 'red', label='The best-fit line')
    plt.scatter(trace_dea['t'], trace_dea['Popen'], c='b', label='The data points')
    plt.legend(loc='best')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(variable_rate + str(variable_rate_val) + '_deaFit.png')
    plt.cla()
    
'''