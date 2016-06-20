# python 3
# pdb parser
# michaladammichalowski@gmail.com
# 25.11.15 - creation
# 10.12.15 - refactor
# EXAMPLE CALL: python3 mm_build_pdbparser.py --old_pdb smth.pdb --new_pdb new_smth.pdb --constrain "name H"

import argparse
import re
import mm_lib_filer as mmfil


class PDBLine:
    def __init__(self, atom_line):
        """
        :param atom_line:
        :return:
        """
        self.line = atom_line
        self.prefix = 'ATOM  '
        self.serial = str(atom_line[6:11])  # shall be int/hex
        self.name = atom_line[11:16]
        self.resName = atom_line[16:20]
        self.chainID = atom_line[20:22]
        self.resSeq = str(atom_line[22:26])  # shall be int/hex
        self.coors = [float(atom_line[30:38]), float(atom_line[38:46]), float(atom_line[46:54])]
        self.occup = float(atom_line[54:60])
        self.temp = float(atom_line[60:66])
        self.extra = atom_line[66:]

    def change_atom(self, temp="default", chain="default"):
        """
        :param temp:
        :param chain:
        :return:
        """
        self.temp = temp if not temp == 'default' else self.temp
        self.chainID = chain if not chain == 'default' else self.chainID

    def write_pdbline(self):
        """
        :param temp:
        :param chain:
        :return:
        """

        # line = '{0:s}{1:5d}{2:s}{3:s}{4:s}{5:4d}{6:s}{7:8.3f}{8:8.3f}{9:8.3f}{10:6.2f}{11:6.2f}{12:s}'.\
        line = '{0:s}{1:s}{2:s}{3:s}{4:s}{5:s}{6:s}{7:8.3f}{8:8.3f}{9:8.3f}{10:6.2f}{11:6.2f}{12:s}'.\
            format\
            (self.prefix,
             self.serial,
             self.name,
             self.resName,
             self.chainID,
             self.resSeq,
             self.line[26:30],
             self.coors[0],
             self.coors[1],
             self.coors[2],
             self.occup,
             self.temp,
             self.extra)

        return line


class PDBFile:
    def __init__(self, old_pdbfile, new_pdbfile):
        """
        :param old_pdbfile:
        :param new_pdbfile:
        :return:
        """
        self.new = new_pdbfile
        with open(old_pdbfile) as infile:
            self.lines = [PDBLine(atom_line) for atom_line in infile.readlines() if atom_line[0:4] == 'ATOM']

    def save(self):
        """
        :return:
        """
        with open(self.new, 'w') as outfile:
            for atom_line in self.lines:
                outfile.write(atom_line.write_pdbline())

    def resegname(self, pattern):
        """
        :param pattern:
        :return:
        """

        # in loop: change mutable objects, not working with immutables
        for atom_line in self.lines:
            atom_line.change_atom(chain=' ' + re.search('{}(.*)'.format(pattern), atom_line.extra.strip()).group(1))

    def remove_name(self, pattern):
        """
        :param pattern:
        :return:
        """

        # in comprehension: avoid removing in loop because of iteration mess
        self.lines = [atom_line for atom_line in self.lines if pattern not in atom_line.name]

    def constrain(self, pattern):
        """
        :param atom_name:
        :return:
        """

        for atom_line in self.lines:
            if pattern in atom_line.extra:
                atom_line.change_atom(temp=1.00)
                print('Constrained: '.format(atom_line.write_pdbline()))
            else:
                atom_line.change_atom(temp=0.00)



parser = argparse.ArgumentParser()
parser.add_argument("-o", "--old_pdb", help="old pdb file extension")
parser.add_argument("-n", "--new_pdb", help="new pdb file extension")
parser.add_argument("-c", "--constrain", help="atom name to constrain")
parser.add_argument("-r", "--remove", help="atom name to remove")
parser.add_argument("-s", "--segname", help="parse segname for chain name, specify prefix to remove")
args = parser.parse_args()


old_files = mmfil.find_files('.', args.old_pdb)
new_files = mmfil.create_outnames(old_files, args.new_pdb)

pdbs = [PDBFile(new, old) for new, old in zip(old_files, new_files)]

for pdb in pdbs:

    if args.constrain:
        pdb.constrain(args.constrain)

    if args.segname:
        pdb.resegname(args.segname)

    if args.remove:
        pdb.remove_name(args.remove)

    pdb.save()
