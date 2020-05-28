# python 3
# concatenates multiple pdb files into one multi model pdb
# ! works fine with modeller 'fit' pdbs
# michaladammichalowski@gmail.com
# EXAMPLE CALL: model_concat.py -m *fit.pdb -o template_type_all_fit.pdb

import argparse


class Models:

    def __init__(self, models_files, output):
        """
        :param models_files: list of pdb file names
        :param output: name of multi model pdb
        """
        self.models_files = models_files
        self.output = output

    def concat(self):
        """
        parses and concatenates pdb files
        """

        models_conc = []

        for model_file in self.models_files:

            header = 'MODEL{:9d}\n'.format(int(model_file[-11:-8]))
            footer = 'ENDMDL\n'
            model_raw = open(model_file, 'r').readlines()
            model = model_raw[1:-1]
            model.insert(0, header)
            model.append(footer)

            models_conc.append(model)

        models_conc_flat = [line for model in models_conc for line in model]
        models_conc_flat.append('END\n')

        with open(self.output, 'w') as f:
            for line in models_conc_flat:
                f.write(line)


parser = argparse.ArgumentParser()
parser.add_argument("-m", "--models_files", dest="models_files", action="store", nargs='+')
parser.add_argument("-o", "--output", dest="output", action="store", type=str)
args = parser.parse_args()

models = Models(args.models_files, args.output)
models.concat()
