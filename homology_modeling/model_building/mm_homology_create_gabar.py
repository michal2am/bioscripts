# python 2
# modeller for GABA execution script
# michaladammichalowski@gmail.com
# ? - creation
# 14.12.15 - refactor
# EXAMPLE CALL: python2 --alnfile a1b2g2_3jad_mm.pir --knows 3JAD_pos --sequence a1b2g2_mm --start_res 8 10 8 10 23 --number 4

from modeller import *
from modeller.automodel import *
from modeller.parallel import *
from mm_homology_gaba_modeller import *
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alnfile", help="alignment file in .pir format")
parser.add_argument("-k", "--knowns", help="known structure name in alignment")
parser.add_argument("-s", "--sequence", help="target sequence name in alignment")
parser.add_argument("-r", "--residues", type=int, nargs="+", help="starting residues of respective chains")
parser.add_argument("-ch", "--chains", nargs="+", help="chain names")
parser.add_argument("-n", "--number", type=int, help="number of models to prepare")
args = parser.parse_args()

j = job()
for core in range(0, 3):
    j.append(local_slave())

# log.verbose()
env = environ(rand_seed=-7777,)                                            # def:-8124, <-2, -50000>
env.io.atom_files_directory = ['.']

gaba_homology = GABAModel(env=env,
                          alnfile=args.alnfile,
                          knowns=args.knowns,
                          sequence=args.sequence,
                          assess_methods=(assess.DOPE, assess.GA341),
                          segments=args.chains,
                          start_res=args.residues)

gaba_homology.starting_model = 1
gaba_homology.ending_model = args.number
gaba_homology.initial_malign3d = False
gaba_homology.final_malign3d = True
gaba_homology.md_level = refine.very_slow
gaba_homology.library_schedule = autosched.slow
gaba_homology.use_parallel_job(j)
gaba_homology.make()
gaba_homology.cluster()
