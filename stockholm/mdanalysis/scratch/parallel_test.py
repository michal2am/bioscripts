# kinda works, but is certainly not fast

import numpy as np
import MDAnalysis
import MDAnalysis.analysis.align
from tqdm import tqdm
import os
import pandas as pd
import argparse
import dask
import dask.multiprocessing
from dask.distributed import Client



parser = argparse.ArgumentParser()
parser.add_argument("--pdb") # full path to pdb
parser.add_argument("--xtc") # xtc file name
parser.add_argument("--dir") # full path to the dir with dir prefix /xxx/yyy/MD
parser.add_argument("--num", type=int) # number of replicas
parser.add_argument("--l_type") # ligand name
parser.add_argument("--l_state") # apo/holo

args = parser.parse_args()

pdb = args.pdb
xtc_f = args.xtc
dir_path = args.dir
dir_total = args.num
ligand_type = args.l_type
ligand_state = args.l_state


pdb_t = '/mnt/as_cephfs/sim/Bicucilline/sim/gromacs/view.pdb'
xtc_t = '/mnt/as_cephfs/sim/Bicucilline/sim/gromacs/MD1/view.xtc'

test_universe = MDAnalysis.Universe(pdb_t, xtc_t)


interface_configs = [['1st_beta/alpha', '((segid A and resid 289) or (segid A and resid 286) or (segid A and resid 265) or (segid B and resid 232))'],
                     ['2nd_beta/alpha', '((segid C and resid 289) or (segid C and resid 286) or (segid C and resid 265) or (segid D and resid 232))'],
                     ['alpha/beta', '((segid B and resid 294) or (segid B and resid 291) or (segid B and resid 270) or (segid C and resid 227))'],
                     ['alpha/gamma', '((segid D and resid 294) or (segid D and resid 291) or (segid D and resid 270) or (segid E and resid 242))'],
                     ['gamma/beta', '((segid E and resid 304) or (segid E and resid 301) or (segid E and resid 280) or (segid A and resid 227))']]


ba1st_atomgroup = test_universe.select_atoms('resname POPC and sphzone 5 ((segid A and resid 289) or (segid A and resid 286) or (segid A and resid 265) or (segid B and resid 232))', updating = True)
ba2nd_atomgroup = test_universe.select_atoms('resname POPC and sphzone 5 ((segid C and resid 289) or (segid C and resid 286) or (segid C and resid 265) or (segid D and resid 232))', updating = True)
ab_atomgroup = test_universe.select_atoms('resname POPC and sphzone 5 ((segid B and resid 294) or (segid B and resid 291) or (segid B and resid 270) or (segid C and resid 227))', updating = True)
ag_atomgroup = test_universe.select_atoms('resname POPC and sphzone 5 ((segid D and resid 294) or (segid D and resid 291) or (segid D and resid 270) or (segid E and resid 242))', updating = True)
gb_atomgroup = test_universe.select_atoms('resname POPC and sphzone 5 ((segid E and resid 304) or (segid E and resid 301) or (segid E and resid 280) or (segid A and resid 227))', updating = True)



def check_presence(atomgroup):
    lipid_atoms_no = len(atomgroup)
    count = 0 if lipid_atoms_no == 0 else 1
    return count


if __name__ == '__main__':

    dask.config.set(scheduler='processes')
    client = Client()
    print(client)

    job_list = []
    for frame in test_universe.trajectory:
        job_list.append(dask.delayed(check_presence(ba1st_atomgroup)))
        job_list.append(dask.delayed(check_presence(ba2nd_atomgroup)))
        job_list.append(dask.delayed(check_presence(ab_atomgroup)))
        job_list.append(dask.delayed(check_presence(ag_atomgroup)))
        job_list.append(dask.delayed(check_presence(gb_atomgroup)))



    result = dask.compute(job_list)
    print(result)
