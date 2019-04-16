import argparse
from mm_homology_sequencereader import Sequence
import pandas as pd


class ScopeEvaluateLocal(Sequence):

    def __init__(self, pir_file, com_seq, sub_seq, sta_res, model_names, chains, energy_files, resi_span, row_names):

        self.models_names = model_names
        self.chains = chains
        self.row_names = row_names

        Sequence.__init__(self, pir_file, com_seq, sub_seq, sta_res)

        self.set_value_csv(com_seq[0], energy_files[1], ['template_ene'])
        self.set_value_csv(com_seq[0], energy_files[0], self.models_names)

        models = self.get_sequence(com_seq[0], num=True, skip=('residue', 'subunit', 'helix', 'strand', 'loop'))

        models_cut = models.loc[(models['position'] > 38) & (models['position'] < 55) & models['subunit'].isin(self.chains)]

        self.both_chains_sums = []
        self.both_chains_singles = []

        for ch, rn in zip(self.chains, self.row_names):

            sep_ch = models_cut.loc[(models['subunit'] == ch)]                                                          # select subunit and create a copy not to mess whole sequence
            singles = sep_ch.copy()
            singles['full_name'] = [rn] * singles.shape[0]

            singles.insert(4, 'mod_avg', singles.loc[:, self.models_names].mean(axis=1))
            print(singles)
            sums = pd.DataFrame(singles.iloc[:, 3:-1].sum()).T                                                          # -1 to drop full name from sum
            sums.index = [rn]

            self.both_chains_sums.append(sums)
            self.both_chains_singles.append(singles)

        self.singles = pd.concat(self.both_chains_singles, axis=0)
        self.singles.to_csv("single.csv")

        self.sums = pd.concat(self.both_chains_sums, axis=0)
        self.sums.index.name = "chain"
        self.sums.to_csv("sumy.csv")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # ScopeEvaluateLocal

    parser.add_argument('-m', '--model_names', nargs='+')
    parser.add_argument('-c', '--chains', nargs='+')
    parser.add_argument('-e', '--energy_files', nargs='+')
    parser.add_argument('-s', '--resi_span', nargs='+')
    parser.add_argument('-n', '--row_names', nargs='+')

    # Sequence
    parser.add_argument('-p',  '--pir_file')
    parser.add_argument('-cs', '--com_seq',   nargs='+')
    parser.add_argument('-ss', '--sub_seq',   nargs='+')
    parser.add_argument('-r',  '--sta_res',   nargs='+', type=int)

    args = parser.parse_args()

    evaluate = ScopeEvaluateLocal(args.pir_file, args.com_seq, args.sub_seq, args.sta_res, args.model_names,
                                  args.chains, args.energy_files, args.resi_span, args.row_names)

