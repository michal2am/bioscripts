# python 3
# pdb parser
# michaladammichalowski@gmail.com
# 25.11.15 - creation
# 10.12.15 - refactor

import argparse


class PDBLine:
    def __init__(self, atom_line):
        """
        :param atom_line:
        :return:
        """
        self.line = atom_line
        self.prefix = 'ATOM  '
        self.serial = int(atom_line[6:11])
        self.name = atom_line[11:16]
        self.resName = atom_line[16:20]
        self.chainID = atom_line[20:22]
        self.resSeq = int(atom_line[22:26])
        self.coors = [float(atom_line[30:38]), float(atom_line[38:46]), float(atom_line[46:54])]
        self.occup = float(atom_line[54:60])
        self.temp = float(atom_line[60:66])
        self.extra = atom_line[66:]

    def write_pdbline(self, temp=''):
        """
        :param temp:
        :return:
        """
        if temp == '':
            temp = self.temp

        line = '{0:s}{1:5d}{2:s}{3:s}{4:s}{5:4d}{6:s}{7:8.3f}{8:8.3f}{9:8.3f}{10:6.2f}{11:6.2f}{12:s}'.format \
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
             temp,
             self.extra)

        return line


class PDBFile:
    def __init__(self, old_pdbfile, new_pdbfile):
        """
        :param old_pdbfile:
        :param new_pdbfile:
        :return:
        """
        self.old = open(old_pdbfile, 'r')
        self.new = open(new_pdbfile, 'w')
        self.old_lines = [PDBLine(atom_line) for atom_line in self.old.readlines() if atom_line[0:4] == 'ATOM']

    def constrain(self, atom_name=''):
        """
        :param atom_name:
        :return:
        """
        for atom_line in self.old_lines:
            if atom_name in atom_line.name:
                self.new.write(atom_line.write_pdbline(temp=1.00))
            else:
                self.new.write(atom_line.write_pdbline())


parser = argparse.ArgumentParser()
parser.add_argument("-o", "--old_pdb", help="old pdb file")
parser.add_argument("-n", "--new_pdb", help="new pdb file")
parser.add_argument("-c", "--constrain", help="atom name to constrain")
args = parser.parse_args()


pdb = PDBFile(args.old_pdb, args.new_pdb)
if args.constrain:
    pdb.constrain(args.constrain)
