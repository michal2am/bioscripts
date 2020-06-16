# python 3
# plots data from data_template_type.csv files
# ! -t specifies template and -v variant order categorization in plots
# michaladammichalowski@gmail.com
# EXAMPLE CALL: python3 model_data_ploters.py -d data* -p plot_name -t temp1 temp2 ... -v var1 var2 ...

import argparse
import pandas as pd
import plotly.express as px


class EvaluatePlots:

    def __init__(self, models_data_files, template_order, variant_order):
        """
        :param models_data_files: standard 'per model' data frame
        """
        read_models = [pd.read_csv(model_data_file) for model_data_file in models_data_files]
        self.models = pd.concat(read_models)
        self.template_order = template_order
        self.variant_order = variant_order

    def plot_distances(self):
        """
        plots two atom pair distances
        """

        fig = px.histogram(self.models, x='distance', color='template',
                           facet_col='variant', facet_row='interface',
                           category_orders={'template': self.template_order, 'variant': self.variant_order},
                           color_discrete_sequence=px.colors.qualitative.G10,
                           height=600, width=600, template='presentation')
        fig.update_layout(barmode='overlay')
        fig.update_traces(opacity=0.66)

        fig.show()
        fig.write_html('histogram_multiModel_distance.html')

    def plot_scores(self):
        """
        plots modellers scores
        """

        fig1 = px.scatter(self.models, x='dope', y='molpdf', hover_name='model',
                          trendline='ols', facet_col='template',
                          color='variant', category_orders={'template': self.template_order,
                                                            'variant': self.variant_order},
                          range_y=[8000, 16000],
                          color_discrete_sequence=px.colors.qualitative.G10,
                          height=500, width=1000, template='presentation',)
        fig1.update_traces(opacity=0.66)
        fig1.show()
        fig1.write_html('scatter_multiModel_score.html')

        fig2 = px.box(self.models, x='template', y='molpdf', hover_name='model',
                      points='all',
                      color='variant', category_orders={'template': self.template_order, 'variant': self.variant_order},
                      range_y=[8000, 16000],
                      color_discrete_sequence=px.colors.qualitative.G10,

                      height=500, width=1000, template='presentation')

        fig2.show()
        fig2.write_html('boxMol_multiModel_score.html')

        fig3 = px.box(self.models, x='template', y='dope', hover_name='model',
                      points='all',
                      color='variant', category_orders={'template': self.template_order, 'variant': self.variant_order},
                      range_y=[-237000, -227000],
                      color_discrete_sequence=px.colors.qualitative.G10,
                      height=500, width=1000, template='presentation')
        fig3.show()
        fig3.write_html('boxDop_multiModel_score.html')


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data", dest="data", action="store", nargs='+')
parser.add_argument("-p", "--plot", dest="plot", action="store", type=str)
parser.add_argument("-t", "--template", dest="template", action="store", nargs='+')
parser.add_argument("-v", "--variant", dest="variant", action="store", nargs='+')

args = parser.parse_args()

evaluater = EvaluatePlots(args.data, args.template, args.variant)

if args.plot == 'score':
    evaluater.plot_scores()
elif args.plot == 'distance':
    evaluater.plot_distances()
