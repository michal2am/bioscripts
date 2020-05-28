# python 3
# parse vmd salt bridge detection files
# ! vmd auto file names should be prefixed with template_variant_ notation
# michaladammichalowski@gmail.com
# EXAMPLE CALL: python model_ss.py -ss *saltbr*

# TODO: model variant hardcoded
# TODO: plot name and properties customization

import argparse
import pandas as pd
import plotly.express as px


class SSBonds:

    def __init__(self, ss_files):
        self.ss_files = ss_files

        self.ss_pairs = self.parse()

    def parse(self):

        ss_pairs = []

        for ss_file in self.ss_files:

            type = ss_file.split('_')[0] + '_' + ss_file.split('_')[1]
            template = ss_file.split('_')[0]
            variant = ss_file.split('_')[1]

            pair_name = ss_file.split('saltbr')[1][1:-4]

            ss_pair = pd.read_csv(ss_file, sep=' ', names=['model', 'distance'])
            ss_pair['model'] = ss_pair['model'] + 1
            ss_pair['type'] = type
            ss_pair['template'] = template
            ss_pair['variant'] = variant
            ss_pair['pair'] = pair_name

            ss_pair = ss_pair.reindex(columns=['model', 'type', 'template', 'variant', 'pair', 'distance'])

            ss_pairs.append(ss_pair)

        ss_pairs = pd.concat(ss_pairs)
        print(ss_pairs)

        averages = pd.pivot_table(ss_pairs, index=['pair'], values='distance', columns='variant')
        averages['doubleCys'] = averages['doubleCys']/averages['wt']
        averages['doubleCysPatch'] = averages['doubleCysPatch']/averages['wt']
        print(averages)

        return ss_pairs

    def plot(self):

        fig = px.box(self.ss_pairs, x='pair', y='distance', hover_name='model',
                      points='all',
                      color='variant', #category_orders={'template': self.template_order, 'variant': self.variant_order},
                      #range_y=[8000, 16000],
                      color_discrete_sequence=px.colors.qualitative.G10,

                      height=750, width=1500, template='presentation')
        fig.show()
        fig.write_html('ss_distances.html')


parser = argparse.ArgumentParser()
parser.add_argument("-ss", "--saltBridges", dest="saltBridges", action="store", type=str, nargs='+')
args = parser.parse_args()

models = SSBonds(args.saltBridges)
models.plot()
