# python 3
# parse rampage critical residues
# michaladammichalowski@gmail.com
# 15.02.16 - creation
#
# EXAMPLE CALL:

import argparse
import re
import collections
import mm_lib_plots as mmplt

class RamachandranResidue:

    def __init__(self, chains, chains_seq, chain, reisd, resname):
        """
        :param chains:
        :param chains_seq:
        :param chain:
        :param reisd:
        :param resname:
        :return:
        """
        self.chains = chains
        self.chains_seq = chains_seq

        self.chain = chain
        self.schain = self.map_subunit()
        self.resid = reisd
        self.resname = resname[1:]  #remove ':' before residue name
        self.name = '{}{}{}'.format(self.resname, self.resid, self.chain)
        self.sname = '{}{}{}'.format(self.resname, self.resid, self.schain)

    def __str__(self):
        """
        :return:
        """
        return '{}{}{}'.format(self.resname, self.resid, self.chain)

    def map_subunit(self):
        """
        :return:
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
        :param model:
        :param rampage_file:
        :param chains:
        :param chains_seq:
        :return:
        """

        self.model = model
        self.file = rampage_file
        self.chains = chains
        self.chains_seq = chains_seq
        self.residues = self.read_residues()

    def read_residues(self):
        """
        :return:
        """
        residues = {'allowed': [], 'outlier': []}
        with open(self.file) as rf:
            rf_residues = [res.rstrip('\n') for res in rf]
            for res in rf_residues:
                if 'Residue' in res:
                    res_info = re.search('\[(.*)\]', res).group(1).split()
                    if 'Allowed' in res:
                        residues['allowed'].append(RamachandranResidue(self.chains, self.chains_seq, *res_info))
                    if 'Outlier' in res:
                        residues['outlier'].append(RamachandranResidue(self.chains, self.chains_seq, *res_info))
        return residues

class RamachandranEvals:

    def __init__(self, evals, frequencies):
        """
        :param evals:
        :return:
        """
        self.evals = evals
        self.freqs = frequencies
        self.bads = {'allowed': self.get_worst('allowed'), 'outlier': self.get_worst('outlier'), \
                     'both': self.get_worst('both')}

    def find_worst(self, which):
        """
        :param which: selected group: allowed or outlier
        :return: dictionary of cumulative amount of which selected residues
        """
        bad_resis = {}
        for model in self.evals:
            for res in model.residues[which]:
                if res.sname in bad_resis.keys():
                    bad_resis[res.sname] += 1
                else:
                    bad_resis.update({res.sname: 1})
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
        :param which:
        :param size:
        :return:
        """
        mmplt.plot_ticker(self.bads[which], 'residue', 'in {} region'.format(which), which, sizex=size)

    def get_print(self, which):
        """
        :param which:
        :param treshold:
        :return:
        """
        counts = collections.Counter(self.bads[which]).most_common()
        print(counts)

parser = argparse.ArgumentParser()
parser.add_argument("--rmf", nargs='+', help="rampage file name")
parser.add_argument("--chains", nargs='+', help="chains in pdb")
parser.add_argument("--chains_seq", nargs='+', help="chains subunit naming")
parser.add_argument("--chains_frq", nargs='+', help="chains frequencies")
args = parser.parse_args()

models = [RamachandranEval('m'+rfile, rfile, args.chains, args.chains_seq) for rfile in args.rmf]
cmodels =RamachandranEvals(models, args.chains_frq)
cmodels.get_plot('outlier', 8.0)
cmodels.get_plot('allowed', 12.0)
cmodels.get_plot('both', 12.0)
cmodels.get_print('outlier')


