# python 3
# parse rampage critical residues
# michaladammichalowski@gmail.com
# 15.02.16 - creation
# 21.03.17 - refactor
#
# EXAMPLE CALL: python3 mm_homology_parse_rama.py  --rmf 31 167 231 242 39  --chains A B C D E --chains_seq B2 A1 B2 A1 Y2 --chains_frq 1 1 1 1 2
import sys
#sys.settrace()
import argparse
import re
import numpy as np
import pandas as pd
import mm_pandas_plot as mmpdplt


class RamachandranResidue:
    def __init__(self, chains_map, chain, resid, resname):
        """
        :param chains_map:
        :param chain: pdb chain
        :param resid: pdb resid
        :param resname: pdb resname
        """
        self.chain, self.schain, self.resid, self.resname = chain, chains_map[chain], resid, resname[1:]
        self.sname = '{}{}{}'.format(self.resname, self.resid, self.schain)


class RamachandranEval:
    def __init__(self, model, rampage_file, chains, chains_seq):
        """
        :param model: model name
        :param rampage_file: rampage outfile name
        :param chains: list of chains as in pdb
        :param chains_seq: list of chain proper naming
        """
        self.model, self.file = model, rampage_file
        self.chains_map = dict(zip(chains, chains_seq))
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
        :return: dict of two lists with residues from allowed and outlier region
        """
        residues = {'allowed': [], 'outlier': []}
        for res in self.read_file('Residue'):
            res_info = re.search('\[(.*)\]', res).group(1).split()
            if 'Allowed' in res:
                residues['allowed'].append(RamachandranResidue(self.chains_map, *res_info))
            if 'Outlier' in res:
                residues['outlier'].append(RamachandranResidue(self.chains_map, *res_info))
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
    def __init__(self, evals, chains_seq, frequencies):
        """
        :param evals:
        :param chains_seq:
        :param frequencies:
        """
        self.freqs = frequencies
        self.chains = chains_seq
        self.evals = evals
        self.bads = pd.DataFrame({'allowed': self.get_worst('allowed'), 'outlier': self.get_worst('outlier')}) \
            .fillna(value=0)

    def get_worst(self, which):
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
                    bad_resis[res.sname] += 1 * equiv
                else:
                    bad_resis.update({res.sname: equiv})
        return bad_resis

    def get_subunit_sum(self):
        """
        :return: total summary for subunits
        """

        # very custom sorting (by subunit and res position)
        self.bads['subunit'] = self.bads.index.map(lambda x: x[-2:])
        self.bads['resnum'] = self.bads.index.map(lambda x: int(''.join([num for num in x[-5:-2] if num.isdigit()])))
        self.bads.sort_values(by=['subunit', 'resnum'], inplace=True, ascending=True)

        # get percentage
        self.bads.loc[:, ['outlier', 'allowed']] = self.bads.loc[:, ['outlier', 'allowed']].apply(
            lambda x: (x / np.max(x)) * 100)

        # greek index
        def greek(resi):
            grk = {'A': r'$\alpha$', 'B': r'$\beta$', 'Y': r'$\gamma$'}[resi[-2]]
            return resi[0:-2] + grk + resi[-1]

        self.bads.index = self.bads.index.map(greek)
        self.bads.index.name = 'residue'

        for subunit in ['A1', 'B2', 'Y2']:

            ploter = mmpdplt.Ploter()
            data = self.bads[self.bads.subunit == subunit].loc[:, ['allowed', 'outlier']]
            print(data.sort_values(by=['outlier', 'allowed'], ascending=False))

            ploter.plot_single('single', data, 'fig_globalrama_' + subunit, [10, 6], y_label='count [%]',
                               kind='bar', rect=(0, 0, 1, 0.85),
                               lines_style={'color': ['grey', 'black'], 'stacked': True},
                               legend_style={'loc': 'best', 'ncol': 2, 'bbox_to_anchor': (0.65, 1.2), 'frameon': True, 'edgecolor': 'black'}
                               )

    def get_model_sum(self):
        """
        :return: total summary of models
        """

        # refactor as fcuk, going pandas
        sums = [[evl.model, evl.summary['favoured'], evl.summary['allowed'], evl.summary['outlier']]
                for evl in self.evals]
        sums = pd.DataFrame(sums, columns=['model', 'favoured', 'allowed', 'outlier']).set_index('model', drop=True)

        # shit done
        sums.sort_values(by=['outlier', 'allowed'], inplace=True, ascending=True)
        print(sums)

print('dupa')
parser = argparse.ArgumentParser()
parser.add_argument("--rmf", nargs='+', help="rampage file name")
parser.add_argument("--chains", nargs='+', help="chains in pdb")
parser.add_argument("--chains_seq", nargs='+', help="chains subunit naming")
parser.add_argument("--chains_frq", nargs='+', type=int, help="chains frequencies")
args = parser.parse_args()

models = [RamachandranEval('m' + rfile, rfile, args.chains, args.chains_seq) for rfile in args.rmf]
cmodels = RamachandranEvals(models, args.chains_seq, args.chains_frq)
cmodels.get_subunit_sum()
cmodels.get_model_sum()
