# Loop refinement of an existing model
from modeller import *
from modeller.automodel import *
from modeller.parallel import *
from mm_homology_gaba_refiner import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--model", help="initial model pdb file")
parser.add_argument("-s", "--sequence", help="sequence name (not relevant)")
parser.add_argument("-n", "--number", type=int, help="number of models to prepare")
parser.add_argument("-l", "--location", help="area to refine (specific residues hardcoded)")
args = parser.parse_args()


j = job()

for core in range(0, 8):
    j.append(local_slave())

log.verbose()
env = environ()
env.io.atom_files_directory = ['.']

m = GABArefine(env, args.model, args.sequence, args.location)

m.loop.starting_model = 1
m.loop.ending_model = args.number
m.loop.md_level = refine.very_slow

m.use_parallel_job(j)
m.make()
