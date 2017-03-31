import argparse
import pandas as pd
import mm_pandas_plot as mmpdplt
import mm_homology_pirreader as mmpirred
import mm_pandas_plot as mmpdplt
from modeller import *
from modeller.scripts import complete_pdb


class EvaluateLocal:

    def __init__(self, name_template, selected_names, template, pir_file, com_seq, sub_seq, sta_res):
        """
        modeller dope evaluation wrapper
        :param name_template: model name prefix
        :param selected_names: model name suffix
        :param template: template name
        """

        self.names = selected_names + ['template']

        self.models = [name_template + selected + '.pdb' for selected in selected_names]
        self.models.append(template)
        self.profiles = ['model_' + selected + '.profile' for selected in selected_names]
        self.profiles.append('template.profile')

        self.pir = mmpirred.Pir(pir_file, com_seq, sub_seq, sta_res)
        self.tempname = com_seq[1]
        self.modname = com_seq[0]

        # self.run()
        self.vals = self.read_profile()
        self.align()
        self.plot()

    def run(self):
        """
        runs modeller dope evaluation and saves results in file 
        """

        log.verbose()  # request verbose output
        env = environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib')  # read topology
        env.libs.parameters.read(file='$(LIB)/par.lib')  # read parameters

        for model, profile in zip(self.models, self.profiles):

            mdl = complete_pdb(env, model)
            s = selection(mdl)
            s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=profile, normalize_profile=True, smoothing_window=15)

    def read_profile(self):
        """
        
        :return: 
        """

        vals = {}

        for name, profile_file in zip(self.names, self.profiles):
            with open(profile_file, 'r') as raw_file:
                vals[name] = []
                for line in raw_file:
                    line = line.strip()
                    if not line.startswith('#') and len(line) > 10:
                        spl = line.split()
                        vals[name].append(float(spl[-1]))
        return vals

    def align(self):
        """
        
        :return: 
        """

        gaps = self.pir.sequences[self.tempname].residue == '-'
        self.pir.sequences[self.tempname].loc[gaps, 'template'] = '-'
        self.pir.sequences[self.tempname].loc[~gaps, 'template'] = self.vals['template']

        gaps = self.pir.sequences[self.modname].residue == '-'
        for comv, val in self.vals.items():
            if comv == 'template': break
            self.pir.sequences[self.modname].loc[gaps, comv] = '-'
            self.pir.sequences[self.modname].loc[~gaps, comv] = val

    def plot(self):
        ploter = mmpdplt.Ploter()
        print( self.pir.sequences['gabaar']['703'])
        ploter.plot_single('single-index', self.pir.sequences['gabaar']['703'], 'fig_localeval', [6, 6], y_label='normalized score',)
                        #   axes_style={'ytvals': [-1, -0.2, -0.5, 0, 0.2, 0.5, 1], 'xtvals': [1000]},
                        #   lines_style={'linestyle': 'None', 'marker': 'o', 'color': ['black', 'grey']},
                        #   legend_style={'loc': 'best', 'ncol': 2}
                        #   )


parser = argparse.ArgumentParser()

# calculate

parser.add_argument("-n", "--nameTemplate")
parser.add_argument("-s", "--selectedNames", nargs='+')
parser.add_argument("-t", "--template")

# pirreader

parser.add_argument('-p',  '--pir_file')
parser.add_argument('-cs', '--com_seq',   nargs='+')
parser.add_argument('-ss', '--sub_seq',   nargs='+')
parser.add_argument('-r',  '--sta_res',   nargs='+', type=int)
parser.add_argument('-sr', '--shift_res', nargs='+', type=int)

args = parser.parse_args()

evaluater = EvaluateLocal(args.nameTemplate, args.selectedNames, args.template, args.pir_file, args.com_seq, args.sub_seq, args.sta_res)
