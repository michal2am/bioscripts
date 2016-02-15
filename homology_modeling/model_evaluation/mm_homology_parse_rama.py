# python 3
# parse rampage critical residues
# michaladammichalowski@gmail.com
# 15.02.16 - creation
#
# EXAMPLE CALL:

import argparse
import mm_lib_plots as mmplt

class RamachandranResidue:

    def __init__(self, chain, reisd, resname):
        self.chain = chain
        self.resid = reisd
        self.resname = resname


class RamachandranEval:

    def __init__(self, model, rampage_file):
        self.model = model
        self.allowed = self.read_residues(rampage_file)['allowed']
        self.outlier = self.read_residues(rampage_file)['outlier']

    def read_residues(self, rampage_file):
        residues = {'allowed': [], 'outlier': []}
        with open(rampage_file) as rf:
            rf_residues = [res.rstrip('\n') for res in rf]
            for res in rf_residues:
                print(res)

        return residues

class RamachandranEvals:
    def __init__(self):
        pass


parser = argparse.ArgumentParser()
parser.add_argument("--rmf", help="rampage file name")
args = parser.parse_args()

test = RamachandranEval('2', args.rmf)
