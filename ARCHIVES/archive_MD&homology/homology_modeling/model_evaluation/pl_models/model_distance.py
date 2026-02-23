# python 3
# calculates distance between two atom pairs in multimodel pdb file
# ! multi model pdb file proper name: template_type_all_fit.pdb, eg. 6i53_doubleCys_all_fit.pdb
# ! parsed pdbs are saved in: data_template_type.csv, eg. data_6i53_doubleCys.csv
# ! all following 'per model' calculations shall be appended to respective data_xxxx_yyyy.csv files as columns
# michaladammichalowski@gmail.com
# EXAMPLE CALL: python3 model_distance.py -m template_type_all_fit.pdb -n model_number -p1 ind1 ind1' -p2 ind2 ind2'

import vmd
import argparse
import numpy as np
import pandas as pd
import plotly.express as px


class MeasureDistances:

    def __init__(self, multi_pdb_file, model_num, pair1, pair2):
        """
        :param multi_pdb_file: pdb file containing all models
        :param model_num: number of models
        :param pair1: indexes of first atom pair
        :param pair2: indexes of second atom pair
        """

        self.multi_pdb_file = multi_pdb_file
        self.model_num = model_num
        self.pair1 = pair1
        self.pair2 = pair2

        self.type = multi_pdb_file[0:-12]
        self.template = self.type.split('_')[0]
        self.variant = self.type.split('_')[1]

        self.distances, self.long_distances = self.parse()

    def parse(self):
        """
        :return: distances between atom pairs in tidy and long format
        """

        molid = vmd.molecule.load('pdb', self.multi_pdb_file)

        dist_dict = {'1st': vmd.measure.bond(atom1=self.pair1[0], atom2=self.pair1[1],
                                             molid=molid, first=0, last=self.model_num - 1),
                     '2nd': vmd.measure.bond(atom1=self.pair2[0], atom2=self.pair2[1],
                                             molid=molid, first=0, last=self.model_num - 1)}
        distance = pd.DataFrame(dist_dict)

        distance['model'] = np.arange(1, self.model_num + 1)
        distance['type'] = self.type
        distance['template'] = self.template
        distance['variant'] = self.variant

        distance = distance.reindex(columns=['model', 'type', 'template', 'variant', '1st', '2nd'])

        long_distance = pd.melt(distance, id_vars=['model', 'type', 'template', 'variant'], value_vars=['1st', '2nd'],
                                var_name='interface', value_name='distance')

        long_distance.to_csv('data_' + self.type + '.csv')
        print(long_distance)

        return distance, long_distance

    def plot_distances(self):
        """
        plots distance histogram
        """

        fig = px.histogram(self.long_distances, x='distance', color='interface',
                           marginal='rug',
                           height=750, width=750, template='presentation')
        fig.update_layout(barmode='overlay')
        fig.update_traces(opacity=0.5)

        # fig.show()
        fig.write_html('box_distance_' + self.type + '.html')


parser = argparse.ArgumentParser()
parser.add_argument("-m", "--multiPdbFile", dest="multiPdbFile", action="store", type=str)
parser.add_argument("-n", "--modelNum", dest="modelNum", action="store", type=int)
parser.add_argument("-p1", "--pair1", dest="pair1", action="store", type=int, nargs='+')
parser.add_argument("-p2", "--pair2", dest="pair2", action="store", type=int, nargs='+')
args = parser.parse_args()

models = MeasureDistances(args.multiPdbFile, args.modelNum, args.pair1, args.pair2)
models.plot_distances()
