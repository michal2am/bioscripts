import pandas as pd
import argparse


class ParseTML:

    def __init__(self, filename):

        # additional info in the file name: TEMPLATE_SUBUNIT_TYPE_extra.txt
        meta_info = filename.split('_')

        # remove garbage from VMD file
        self.parsed = pd.read_csv(filename, sep=' ', header=None)
        self.parsed.drop(axis=0, index=range(9), inplace=True)
        self.parsed.drop(axis=1, columns=2, inplace=True)
        self.parsed.rename(columns={0: 'position', 1: 'subunit', 3: 'model', 4: 'rmsf'}, inplace=True)

        # cast correct types
        self.parsed['position'] = self.parsed.position.astype(int)
        self.parsed['model'] = self.parsed.model.astype(int)
        self.parsed['rmsf'] = self.parsed.rmsf.astype(float)

        # add info from file name
        self.parsed['template'] = meta_info[0]
        self.parsed['subunit'] = meta_info[1]
        self.parsed['type'] = meta_info[2]

        # tidy up column order
        new_order = ['position', 'template', 'subunit', 'type', 'model', 'rmsf']
        self.parsed = self.parsed.reindex(columns=new_order)

        print(self.parsed)

        self.parsed.to_csv(filename.strip('.txt') + '_parsed.csv')

        self.parsed_heatmap = self.parsed.pivot(index='position', columns='model', values='rmsf')

        print(self.parsed_heatmap)


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file_name')
args = parser.parse_args()


parse = ParseTML(args.file_name)

