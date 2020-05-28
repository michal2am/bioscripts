# python 3
# parses and plots various modeller scores values from log files
# ! modeller log proper name: template_type.log, eg. 6i53_doubleCys.log
# ! parsed logs are saved in: data_template_type.csv, eg. data_6i53_doubleCys.csv
# ! all following 'per model' calculations shall be appended to respective data_xxxx_yyyy.csv files as columns
# michaladammichalowski@gmail.com
# EXAMPLE CALL: python3 mm_homology_eval_global.py -l template_type.log -n model_number -t float_percent_number

import subprocess as sp
import argparse
import numpy as np
import pandas as pd
import plotly.express as px


class EvaluateGlobal:

    def __init__(self, log_file, model_num, threshold):
        """
        :param log_file: modeller log file
        :param model_num: number of models
        :param threshold: percent of top models to select
        """

        self.log_file, self.model_num, self.thresh = log_file, model_num, int(threshold * model_num)
        self.type = log_file.split('.')[0]
        self.all_models, self.best_molpdf, self.best_dope, self.best_common, self.normals = self.parse()

    def parse(self):
        """
        :return: all models, best x% molpdf, best x% dope, best common, best common normalized
        """

        model_list = [model.split() for model in sp.check_output(['grep', 'Summary of successfully produced models:',
                      self.log_file, '-A', str(self.model_num+2)])
                      .decode('utf-8').split('\n')[3:-1]]
        model_list = [[int(model[0][17:20]), float(model[1]), float(model[2])] for model in model_list]

        all_models = pd.DataFrame(model_list, columns=['model', 'molpdf', 'dope'])

        best_molpdf = all_models.sort_values(by=['molpdf']).iloc[0:self.thresh]
        best_dope = all_models.sort_values(by=['dope']).iloc[0:self.thresh]

        common = best_molpdf.index.intersection(best_dope.index)
        c_by_molpdf = all_models.loc[common].sort_values(by=['molpdf'])
        n_by_model = all_models.loc[:, ['molpdf', 'dope']].apply(lambda x: (x - np.mean(x))/(np.min(x) - np.max(x)))

        all_models['type'] = self.type
        all_models['template'] = self.type.split('_')[0]
        all_models['variant'] = self.type.split('_')[1]

        all_models = all_models.reindex(columns=['model', 'type', 'template', 'variant', 'molpdf', 'dope'])
        all_models.to_csv('data_' + self.type + '.csv')

        print('Top by molpdf:\n', best_molpdf)
        print('Top by dope:\n', best_dope)
        print('Top common:\n', c_by_molpdf)

        return all_models, best_molpdf, best_dope, c_by_molpdf, n_by_model

    def plot_scores(self):
        """
        plots dope vs molpdf scatter
        """

        fig = px.scatter(self.all_models, x='dope', y='molpdf', hover_name='model',
                         trendline='ols', marginal_x='rug', marginal_y='rug',
                         title=self.type,
                         height=750, width=750, template='presentation')
        fig.show()
        fig.write_html('scatter_modelScore_' + self.type + '.html')


parser = argparse.ArgumentParser()
parser.add_argument("-l", "--logFile", dest="logFile", action="store", type=str)
parser.add_argument("-n", "--modelNum", dest="modelNum", action="store", type=int)
parser.add_argument("-t", "--threshold", dest="threshold", action="store", type=float)
args = parser.parse_args()

models = EvaluateGlobal(args.logFile, args.modelNum, args.threshold)
models.plot_scores()
