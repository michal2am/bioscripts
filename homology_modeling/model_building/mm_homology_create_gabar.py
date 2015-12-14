# python 2
# modeller for GABA execution script
# michaladammichalowski@gmail.com
# ? - creation
# 14.12.15 - refactor
# EXAMPLE CALL: python2 --alnfile a1b2g2_3jad_mm.pir --knows 3JAD_pos --sequence a1b2g2_mm --start_res 8 10 8 10 23 --number 4

from modeller import *
from modeller.parallel import *
from mm_homology_gaba_modeller import *
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alnfile", help="alignment file in .pir format")
parser.add_argument("-k", "--knows", help="known structure name in alignment")
parser.add_argument("-s", "--sequence", help="target sequence name in alignment")
parser.add_argument("-r", "--residues", nargs="+", help="starting residues of respective chains")
parser.add_argument("-n", "--number", help="number of models to prepare")
args = parser.parse_args()


j = job()
j.append(local_slave())
j.append(local_slave())
j.append(local_slave())
j.append(local_slave())

log.verbose()   
env = environ() 
env.io.atom_files_directory = ['.']

a = GABAModel(env,
              alnfile=args.alnfile,
              knowns=args.knows,
              sequence=args.sequence,
              assess_methods=(assess.DOPE, assess.GA341),
              segments=['A', 'B', 'C', 'D', 'E'],
              start_res=args.residues)

a.starting_model = 1
a.ending_model = args.number
a.initial_malign3d = False
a.final_malign3d = True
a.md_level = refine.very_slow
a.library_schedule = autosched.slow  
a.use_parallel_job(j)
a.cluster()
a.make()                            

