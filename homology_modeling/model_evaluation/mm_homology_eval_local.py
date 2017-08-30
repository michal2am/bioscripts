import argparse
import mm_pandas_plot as mmpdplt
from modeller import *
from modeller.scripts import complete_pdb
from mm_homology_sequencereader import Sequence
import pandas as pd


class EvaluateLocal(Sequence):

    def __init__(self, name_template, selected_names, template, pir_file, com_seq, sub_seq, sta_res, run):
        """
        modeller dope evaluation wrapper
        :param name_template: model name prefix
        :param selected_names: model name suffix
        :param template: template file name
        :param pir_file: pir file name
        :param com_seq: name of complete sequence
        :param sub_seq: name of subunitp
        :param sta_res: starting residue numbers
        """

        Sequence.__init__(self, pir_file, com_seq, sub_seq, sta_res)

        self.set_value_csv('gabaar', 'csv_loopnaming_gabaarlinkershort.csv', ['strand', 'helix', 'loop'])

        self.temp_name = self.com_seq[1]
        self.mod_name = self.com_seq[0]

        self.names = selected_names + ['template']
        self.models = [name_template + selected + '.pdb' for selected in selected_names]
        self.models.append(template)
        self.profiles = ['model_' + selected + '.profile' for selected in selected_names]
        self.profiles.append('template.profile')

        if run: self.run()
        self.read_profile()
        self.plot()

    def run(self):
        """
        runs modeller dope evaluation and saves results in file 
        """

        log.verbose()
        env = environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')

        for model, profile in zip(self.models, self.profiles):

            mdl = complete_pdb(env, model)
            s = selection(mdl)
            s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=profile, normalize_profile=True, smoothing_window=15)

    def read_profile(self):
        """ 
        reads modeller dope evaluation files and adds them to Sequence
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

        temp_val = {'template': vals.pop('template')}
        Sequence.set_value(self, self.temp_name, temp_val)
        Sequence.set_value(self, self.mod_name, vals)



    def plot(self):
        """
        plots results for respective subunits of template and models
        """

        ploter = mmpdplt.Ploter()

        self.get_sequence(self.mod_name).to_csv('csv_localeval_' + self.mod_name + '.csv')
        self.get_sequence(self.temp_name).to_csv('csv_localeval_' + self.temp_name + '.csv')

        for seq_name in zip(self.sub_seq[0:self.sub_num], self.sub_seq[self.sub_num:]):                                 # model & template subunit pairs

            to_plot = [self.get_sequence(self.temp_name, sub_name=seq_name[1], num=True, skip=('residue', 'position', 'subunit', 'strand', 'helix', 'loop')),
                       self.get_sequence(self.mod_name, sub_name=seq_name[0], num=True, skip=('residue', 'position', 'subunit', 'strand', 'helix', 'loop'))]
            print(to_plot)

            title = {'a': r'$\alpha$', 'b': r'$\beta$', 'g': r'$\gamma$'}[seq_name[0][0]] + seq_name[0][1:]

            ploter.plot_single('multi', [frame.drop('position', axis=1) for frame in to_plot], 'fig_localeval_{}'.format(seq_name[0]), [10, 6],
                               y_label='normalized score', x_label='residue',
                               axes_style={'xtlabs':  to_plot[1]['position'][1::50],
                                           'xtvals': to_plot[1].index.values[1::50]},
                               lines_style={'linestyle': '-', 'marker': 'o', 'title': title},
                               legend_style={'bbox_to_anchor':(0.5, .5), 'ncol': 3, 'frameon': True, 'edgecolor': 'black'},
                               annotation={'features': ['strand', 'helix', 'loop']}
                               )


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # EvaluateLocal

    parser.add_argument("-n", "--nameTemplate")
    parser.add_argument("-s", "--selectedNames", nargs='+')
    parser.add_argument("-t", "--template")

    # Sequence

    parser.add_argument('-p',  '--pir_file')
    parser.add_argument('-cs', '--com_seq',   nargs='+')
    parser.add_argument('-ss', '--sub_seq',   nargs='+')
    parser.add_argument('-r',  '--sta_res',   nargs='+', type=int)
    parser.add_argument('-sr', '--shift_res', nargs='+', type=int)
    parser.add_argument("-rn", "--run", dest='run', action='store_true')
    parser.add_argument("-nrn", "--no_run", dest='run', action='store_false')

    parser.set_defaults(run=False)

    args = parser.parse_args()

    evaluater = EvaluateLocal(args.nameTemplate, args.selectedNames, args.template, args.pir_file, args.com_seq,
                              args.sub_seq, args.sta_res, args.run)
