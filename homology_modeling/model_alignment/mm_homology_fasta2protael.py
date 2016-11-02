import re


class Alignment:

    def __init__(self, f):
        self.file = f
        self.sequences = self.read_aln()

    def read_aln(self):
        with open(self.file) as aln_file:
            alns = aln_file.readlines()
            n_alns = []
            s_alns = []
            for line in alns:
                if line[0] == '>':
                    n_alns.append(re.search('>(.*)/.*', line).group(1))
                    sequence = ""
                else:
                    sequence += line.replace('\n', '').replace('\r', '')
        return p_alns

    def write_prot(self):
        for name in self.sequences:
            print(name)
            print(Protael.scheme(name, self.sequences[name], 'none'))


class Protael:

    @staticmethod
    def scheme(name, sequence, pdb):
        return """""{
                    "CS": "Clustal",
                    "id": "{0}",
                    "label": "{0}",
                    "sequence": "{1}",
                    "start": 6,
                    "properties": {
                                  "pdbid": "{2}",
                                  "PDB structure": "{2}"
                                  }
                    },""".format(name, sequence, pdb)


test_aln = Alignment('corrected_mid_reorder_nogaps_all_result.fasta')
test_aln.write_prot()