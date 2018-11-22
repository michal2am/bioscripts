import pandas as pd
import argparse as ap
import copy as cp


class Sequence:

    def __init__(self, pir_file, com_seq, sub_seq, sta_res):
        """
        stores pir type alignment with additional data
        :param pir_file: pir file name
        :param com_seq: name of complete sequence
        :param sub_seq: name of subunit
        :param sta_res: starting residue numbers
        """

        self.pir_file = pir_file
        self.com_seq, self.sub_seq, self.sta_res = com_seq, sub_seq, sta_res
        self.sub_num = int(len(self.sub_seq)/len(self.com_seq))  #WTF? skąd to dzielenie?
        #self.sub_num = int(len(self.sub_seq))

        """
        main attribute:
        dictionary of complete sequences
        dictionary key is a name of complete sequence
        dictionary item is a dataframe, basic columns:
            residue number (index)
            columns: residue (one letter code), subunit, position (physiological index)
        """
        self.sequences = self.read()

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

        sequences = (dict(zip(self.sub_seq, sequences)))                                                                # adds 'subunit' column
        for name, seq in sequences.items():
            seq['subunit'] = name

        def numbering(sequences):
            """
            assigns physiological residue numbering to 'position' column
            """

            for seq, sta in zip(sequences.items(), self.sta_res):
                gaps = seq[1].residue == '-'
                seq[1].loc[gaps, 'position'] = '-'
                seq[1].loc[~gaps, 'position'] = pd.np.arange(sta, (~gaps).sum() + sta)

        def concatenating(self, sequences):
            """
            concatenates subunits of same complete sequence
            :return: dictionary of concatenated subunit sequences
            """

            sub_cat = {com: self.sub_seq[i * self.sub_num:i * self.sub_num + self.sub_num] for i, com in
                       enumerate(self.com_seq)}
            print(sub_cat)
            sequences = {com: pd.concat([sequences[sub] for sub in sub_cat[com]], ignore_index=True)
                         for com in sub_cat}

            for seq in sequences: sequences[seq].index.name = 'residue number'

            return sequences

        numbering(sequences)
        sequences = concatenating(self, sequences)

        return sequences

    def numerize(self, skip=('residue', 'position', 'subunit')):
        """
        convert values in sequences to numeric
        :param skip: column names to skip (list)
        :return: numerized copy of sequences
        """

        cp_sequences = cp.deepcopy(self.sequences)

        for sequence, values in cp_sequences.items():
            for column in values.ix[:, values.columns.difference(skip)]:
                values[column] = pd.to_numeric(values[column], errors='coerce')

        return cp_sequences

    def get_sequence(self, seq_name, sub_name=False, num=False, skip=('residue', 'position', 'subunit')):
        """
        select subunit from all sequence dictionary
        :param seq_name: name of complete sequence (dict key)
        :param sub_name: name of subunit (df column name)
        :param num: numerize columns (bool)
        :param skip: columns to skip at numerizing (list)
        :return: dataframe of selected subunit only
        """
        if sub_name:
            sub_mask = self.sequences[seq_name]['subunit'] == sub_name
        else:
            sub_mask = self.sequences[seq_name].index.values

        if num:
            return self.numerize(skip)[seq_name].loc[sub_mask]
        else:
            return self.sequences[seq_name].loc[sub_mask]

    def set_value(self, seq_name, values):
        """
        sets additional columns for selected complete sequence
        :param seq_name: name of complete sequence (dict key)
        :param values: values to add (dict)
        """
        print(seq_name)
        print(values)
        gaps = self.sequences[seq_name].residue == '-'
        #print(gaps)
        print(self.sequences[seq_name])

        for col_name, col_val in values.items():
            self.sequences[seq_name].loc[gaps, col_name] = '-'
            self.sequences[seq_name].loc[~gaps, col_name] = col_val

    def set_value_csv(self, seq_name, csv, headers):
        """
        sets additional columns for selected complete sequence from csv file
        :param seq_name: ame of complete sequence (dict key)
        :param csv: csv file name
        :param headers: headers in csv to append
        :return: 
        """
        datafile = pd.read_csv(csv)
        self.set_value(seq_name, {column: datafile[column].values for column in headers})



    def shift_numbering(self):
        pass

if __name__ == '__main__':

    parser = ap.ArgumentParser()
    parser.add_argument('-p',  '--pir_file')
    parser.add_argument('-cs', '--com_seq',   nargs='+')
    parser.add_argument('-ss', '--sub_seq',   nargs='+')
    parser.add_argument('-r',  '--sta_res',   nargs='+', type=int)
    parser.add_argument('-sr', '--shift_res', nargs='+', type=int)
    parser.add_argument('-s',  '--shift',     nargs='+', type=int)

    args = parser.parse_args()
    pir = Sequence(args.pir_file, args.com_seq, args.sub_seq, args.sta_res)
    pir.get_sequence('gabaar').to_csv('test.csv')








