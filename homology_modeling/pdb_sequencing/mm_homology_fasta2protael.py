import re
import json
from collections import OrderedDict


class Alignment:

    def __init__(self, f):

        self.file = f
        self.sequences = self.read_aln()

    def read_aln(self):

        with open(self.file) as aln_file:
            alns = aln_file.readlines()
            sequences = OrderedDict()
            for line in alns:
                if line[0] == '>':
                    name = (re.search('>(.*)/.*', line).group(1))
                    sequences[name] = ''
                else:
                    sequences[name] += line.replace('\n', '').replace('\r', '')
        return sequences

    def write_prot(self):

        json_str = []

        def l2grekk(n, id):
            codes = {'A': '\u03B1', 'B': '\u03B2', 'G': '\u03B3', 'D': '\u03B4', 'R': '\u03C1', 'E': '\u03B5',
                    'T': '\u03B8', 'Pi': '\u03C0'}

            n_n = n.split('_')
            try:
                for code in codes:
                    n_n[id] = n_n[id].replace(code, codes[code])
                n_n = ' '.join(n_n)
            except IndexError:
                n_n = ' '.join(n_n)

            return n_n

        for name in self.sequences:

            try:
                pdb = re.search('[a-z, A-Z]*_([A-Z, 0-9]{4}).*', name).group(1)
                n_name = l2grekk(name, 2)

            except AttributeError:
                pdb = 'none'
                n_name = l2grekk(name, 1)

            print(name, n_name)

            full = OrderedDict((('CS', 'Clustal'), ('id', n_name), ('label', n_name), ('sequence', self.sequences[name]),
                               ('start', 6), ('properties', OrderedDict((('pdbid', pdb), ('PDB structure', pdb))))))
            json_str.append(full)

            with open('tmp', 'w') as f:
                json.dump(json_str, f, indent=5)

test_aln = Alignment('corrected_mid_reorder_nogaps_all_result.fasta')
test_aln.write_prot()
