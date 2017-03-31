import pandas as pd
# import argparse as ap


class Pir:

    def __init__(self, pir_file, com_seq, sub_seq, sta_res):

        self.pir_file = pir_file
        self.com_seq, self.sub_seq, self.sta_res = com_seq, sub_seq, sta_res

        self.sub_sequences = self.read()
        self.numbering()
        self.sequences = self.concatenating()

    def read(self):
        """
        reads pir file
        :return: dictionary of read molecules' sequences in data frames in 'residue' column
        """

        sequences = []

        with open(self.pir_file, 'r') as raw_file:
            for line in raw_file:
                line = line.strip()
                if not any(starter in line for starter in ('>', 'sequence', 'structure')):
                    sequences.append(line)

        sequences = [seq.split('/') for seq in ''.join(sequences).split('*')[0:-1]]                                     # join & split by molecule and subunits
        sequences = [pd.DataFrame(list(seq), columns=['residue']) for sub_seq in sequences for seq in sub_seq]          # flat and string to df

        sequences = (dict(zip(self.sub_seq, sequences)))

        return sequences

    def numbering(self):
        """
        assigns physiological residue numbering to 'position' column
        """

        for seq, sta in zip(self.sub_sequences.items(), self.sta_res):

            gaps = seq[1].residue == '-'
            seq[1].loc[gaps, 'position'] = '-'
            seq[1].loc[~gaps, 'position'] = pd.np.arange(sta, (~gaps).sum() + sta)

    def concatenating(self):
        """
        concatenates subunits of same protein
        :return: dictionary of concatenated subunit sequences
        """

        sub_num = int(len(self.sub_seq)/len(self.com_seq))

        sub_cat = {com: self.sub_seq[i*sub_num:i*sub_num+sub_num] for i, com in enumerate(self.com_seq)}
        sequences = {com: pd.concat([self.sub_sequences[sub] for sub in sub_cat[com]], ignore_index=True)
                     for com in sub_cat}

        return sequences

    def shift_numbering(self):
        pass

'''

parser = ap.ArgumentParser()
parser.add_argument('-p',  '--pir_file')
parser.add_argument('-cs', '--com_seq',   nargs='+')
parser.add_argument('-ss', '--sub_seq',   nargs='+')
parser.add_argument('-r',  '--sta_res',   nargs='+', type=int)
parser.add_argument('-sr', '--shift_res', nargs='+', type=int)
parser.add_argument('-s',  '--shift',     nargs='+', type=int)

args = parser.parse_args()
pir = Pir(args.pir_file, args.com_seq, args.sub_seq, args.sta_res)

'''







