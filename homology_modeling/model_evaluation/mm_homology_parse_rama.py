# python 3
# parse rampage critical residues
# michaladammichalowski@gmail.com
# 15.02.16 - creation
#
# EXAMPLE CALL: python3 mm_homology_parse_rama.py --model glucl_3JADcl  --rmf 31 167 231 242 39  --chains A B C D E --chains_seq B2 A1 B2 A1 Y2 --chains_frq 1 1 1 1 2


import argparse
import re
import collections
import mm_lib_plots as mmplt
import mm_lib_analysis as mmanl


class RamachandranResidue:

    def __init__(self, chains, chains_seq, chain, reisd, resname):
        """
        :param chains: list of chains as in pdb
        :param chains_seq: list of chain proper naming
        :param chain: pdb chain
        :param reisd: pdb resid
        :param resname: pdb resname
        :return:
        """
        self.chains = chains
        self.chains_seq = chains_seq

        self.chain = chain
        self.schain = self.map_subunit()
        self.resid = reisd
        self.resname = resname[1:]  # remove ':' before residue name
        self.name = '{}{}{}'.format(self.resname, self.resid, self.chain)
        self.sname = '{}{}{}'.format(self.resname, self.resid, self.schain)

    def __str__(self):
        """
        :return: resnameResidChain in proper naming manner
        """
        return '{}{}{}'.format(self.resname, self.resid, self.schain)

    def map_subunit(self):
        """
        :return: maps pdb chains to proper chain naming
        """
        subunit = ''
        for ch, chs in zip(self.chains, self.chains_seq):
            if re.search(r"\b{0}\b".format(ch), self.chain):
                subunit = chs
                break
        return subunit


class RamachandranEval:

    def __init__(self, model, rampage_file, chains, chains_seq):
        """
        :param model: model name
        :param rampage_file: rampage outfile name
        :param chains: list of chains as in pdb
        :param chains_seq: list of chain proper naming
        :return:
        """

        self.model = model
        self.file = rampage_file
        self.chains = chains
        self.chains_seq = chains_seq
        self.residues = self.read_residues()
        self.summary = self.read_summary()

    def read_file(self, phrase):
        """
        :param phrase: phrase to look for in file
        :return: list of file lines containing phrase
        """
        parsed = []
        with open(self.file) as rf:
            rf_residues = [res.rstrip('\n') for res in rf]
            for res in rf_residues:
                if phrase in res:
                    parsed.append(res)
        return parsed

    def read_residues(self):
        """
        :return: dict of two lists with resiudes from allowed and outlier region
        """
        residues = {'allowed': [], 'outlier': []}
        for res in self.read_file('Residue'):
            res_info = re.search('\[(.*)\]', res).group(1).split()
            if 'Allowed' in res:
                residues['allowed'].append(RamachandranResidue(self.chains, self.chains_seq, *res_info))
            if 'Outlier' in res:
                residues['outlier'].append(RamachandranResidue(self.chains, self.chains_seq, *res_info))
        return residues

    def read_summary(self):
        """
        :return: dict fot total percent of favoured, allowed and outlier residues
        """
        summary = {'favoured': '', 'allowed': '', 'outlier': ''}
        for res in self.read_file('Number of residues in'):
            persc = re.search(':.*\((.*%)\)', res).group(1)
            if 'favoured' in res:
                summary['favoured'] = persc
            if 'allowed' in res:
                summary['allowed'] = persc
            if 'outlier' in res:
                summary['outlier'] = persc
        return summary


class RamachandranEvals:

    def __init__(self, model, evals, chains_seq, frequencies):
        """
        :param model:
        :param evals:
        :param chains_seq:
        :param frequencies:
        :return:
        """
        self.gen_model = model
        self.evals = evals
        self.freqs = frequencies
        self.chains = chains_seq
        self.bads = {'allowed': self.get_worst('allowed'), 'outlier': self.get_worst('outlier'),
                     'both': self.get_worst('both')}

    def find_worst(self, which):
        """
        :param which: selected group: allowed or outlier
        :return: dictionary of cumulative amount of which selected residues
        """
        bad_resis = {}
        for model in self.evals:
            for res in model.residues[which]:
                chain_ind = self.chains.index(res.schain)
                equiv = self.freqs[chain_ind]
                if res.sname in bad_resis.keys():
                    bad_resis[res.sname] += 1*equiv
                else:
                    bad_resis.update({res.sname: equiv})
        return bad_resis

    def get_worst(self, which):
        """
        find_worst wrapper
        :param which: selected group: allowed, outlier or both
        :return: dictionary of cumulative amount of which selected residues (wrapper for find_worst)
        """
        if which in ['allowed', 'outlier']:
            return self.find_worst(which)
        if which == 'both':
            allowed = self.find_worst('allowed')
            outlier = self.find_worst('outlier')
            both = allowed.copy()
            both.update(outlier)
            return both

    def get_plot(self, which, size):
        """
        :param which: allowed, outlier or both residues
        :param size: plot of number of residues in selected region
        :return:
        """
        mmplt.plot_ticker(self.bads[which], 'residue', 'in {} region'.format(which), which, sizex=size)

    def get_print_worst(self, which):
        """
        :param which: allowed, outlier or both residues
        :return: numerical summary of number of residues in selected region
        """
        print('Residues in {} region:'.format(which))
        for res in collections.Counter(self.bads[which]).most_common():
            print(*res, sep=' ', end=', ')
        print('\n')

    def get_print_sum(self):
        """
        :return: total summary of models
        """
        sums = [[evl.model, evl.summary['favoured'], evl.summary['allowed'], evl.summary['outlier']]
                for evl in self.evals]
        sums = mmanl.get_sort(sums, 3)
        cols = '   model favoured allowed outlier'
        sep = '-' * len(cols)
        print('{}\n{} rampage evaluation\n{}\n{}\n{}'.format(sep, self.gen_model, sep, cols, sep))
        for sm in sums:
            print('{}  {}\t{}\t{}\t{}'.format(*sm))
        print(sep)

parser = argparse.ArgumentParser()
parser.add_argument("--model", help="general name of the model")
parser.add_argument("--rmf", nargs='+', help="rampage file name")
parser.add_argument("--chains", nargs='+', help="chains in pdb")
parser.add_argument("--chains_seq", nargs='+', help="chains subunit naming")
parser.add_argument("--chains_frq", nargs='+', type=int, help="chains frequencies")
args = parser.parse_args()

models = [RamachandranEval('m'+rfile, rfile, args.chains, args.chains_seq) for rfile in args.rmf]
cmodels = RamachandranEvals(args.model, models, args.chains_seq, args.chains_frq)
cmodels.get_plot('outlier', 5.0)
cmodels.get_plot('allowed', 12.0)
cmodels.get_plot('both', 12.0)
cmodels.get_print_worst('allowed')
cmodels.get_print_worst('outlier')
cmodels.get_print_worst('both')
cmodels.get_print_sum()
