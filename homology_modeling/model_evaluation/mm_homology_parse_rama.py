# python 3
# parse rampage critical residues
# michaladammichalowski@gmail.com
# 15.02.16 - creation
#
# EXAMPLE CALL:

import argparse
import re
import mm_lib_plots as mmplt

class RamachandranResidue:

    def __init__(self, chain, reisd, resname):
        """
        :param chain:
        :param reisd:
        :param resname:
        :return:
        """

        self.chain = chain
        self.resid = reisd
        self.resname = resname[1:]  #remove ':' before residue name
        self.name = '{}{}{}'.format(self.resname, self.resid, self.chain)

    def __str__(self):
        """
        :return:
        """

        return '{}{}{}'.format(self.resname, self.resid, self.chain)


class RamachandranEval:

    def __init__(self, model, rampage_file):
        """
        :param model:
        :param rampage_file:
        :return:
        """

        self.model = model
        self.file = rampage_file
        self.residues = self.read_residues()

    def read_residues(self):
        """
        :param rampage_file:
        :return:
        """

        residues = {'allowed': [], 'outlier': []}
        with open(self.file) as rf:
            rf_residues = [res.rstrip('\n') for res in rf]
            for res in rf_residues:
                if 'Residue' in res:
                    res_info = re.search('\[(.*)\]', res).group(1).split()
                    if 'Allowed' in res:
                        residues['allowed'].append(RamachandranResidue(*res_info))
                    if 'Outlier' in res:
                        residues['outlier'].append(RamachandranResidue(*res_info))
        return residues

class RamachandranEvals:

    def __init__(self, evals):
        self.evals = evals

    def find_worst(self):
        bad_resis = {}
        for model in self.evals:
            for res in model.residues['allowed']:
                if res.name in bad_resis.keys():
                    bad_resis[res.name] += 1
                else:
                    bad_resis.update({res.name: 1})



parser = argparse.ArgumentParser()
parser.add_argument("--rmf", help="rampage file name")
args = parser.parse_args()

test = RamachandranEval('31', args.rmf)
for res in test.allowed: print(res)