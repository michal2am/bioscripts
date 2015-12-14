# python 2
# script for transcripting pdb to fasta
# michaladammichalowski@gmail.com
# ? - creation
# 14.12.15 - refactor
# EXAMPLE CALL: mm_pdb2fasta.py --pdbFile smth.pdb --pdbCode smth --seqName smth_seq

import argparse
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder


def getSeq(structName, chain, ch, seqFile):
    ppb = PPBuilder()
    polypeptide = ppb.build_peptides(chain)

    resNumbers = []
    for residue in chain:
        if residue.id[0] == ' ':
            resNumbers.append(residue.id[1])

    SeqSta = resNumbers[0]
    SeqSto = resNumbers[-1]
    PdbSeq = ''

    for peptide in polypeptide:
        PdbSeq += str(peptide.get_sequence())
    print '>{0} chain-{1}'.format(structName, ch)
    print PdbSeq

    seqFile.write('>{0} chain-{1}\n'.format(structName, ch))
    seqFile.write(PdbSeq + '\n')


def parsePdb(pdbFile, pdbCode, seqName):
    structure = PDBParser().get_structure(pdbCode, pdbFile)
    chains = ['A', 'B', 'C', 'D', 'E']
    seqFile = open(seqName + '.fasta', 'w')

    for ch in chains:
        chain = structure[0][ch]
        getSeq(pdbCode, chain, ch, seqFile)

    seqFile.close()


parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pdbFile", dest="pdbFile", action="store")
parser.add_argument("-c", "--pdbCode", dest="pdbCode", action="store")
parser.add_argument("-n", "--seqName", dest="seqName", action="store")
args = parser.parse_args()

parsePdb(args.pdbFile, args.pdbCode, args.seqName)
