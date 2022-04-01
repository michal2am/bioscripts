import numpy as np
import MDAnalysis
import MDAnalysis.analysis.align
from tqdm import tqdm
import os
import pandas as pd

pdb = 'step5_input_chains.pdb'
xtc_f = 'all_runsteps.xtc'
dir_name = 'sys'
ligand_type = 'etomidate'
ligand_state = 'apo'

all_complete = pd.DataFrame()

for i in range(1, 5):

    print('Going to sys{}'.format(i))
    xtc = dir_name + str(int(i)) + '/gromacs/' + xtc_f

    if os.path.isfile(xtc):

        U = MDAnalysis.Universe(pdb, xtc)

        first_bs = '(segid A and resid 289) or (segid A and resid 286) or (segid A and resid 265) ' \
                   'or (segid B and resid 232) or (segid B and resid 233)'

        first_lipid = U.select_atoms('((resname POPC) and sphzone 5.0 {})'.format(first_bs), updating=True)

        second_bs = '(segid C and resid 289) or (segid C and resid 286) or (segid C and resid 265) ' \
                    'or (segid D and resid 232) or (segid D and resid 233)'

        second_lipid = U.select_atoms('((resname POPC) and sphzone 5.0 {})'.format(second_bs), updating=True)

        first_occupied = []
        second_occupied = []

        for ts in tqdm(U.trajectory):
            first_occupied.append(0) if len(first_lipid) == 0 else first_occupied.append(1)
            second_occupied.append(0) if len(second_lipid) == 0 else second_occupied.append(1)

        first_complete = pd.DataFrame()

        first_complete['occupied'] = first_occupied
        first_complete['system'] = dir_name + str(int(i))
        first_complete['interface'] = '1st'
        first_complete['ligand_type'] = ligand_type
        first_complete['ligand_state'] = ligand_state

        all_complete = pd.concat([all_complete, first_complete])

        second_complete = pd.DataFrame()

        second_complete['occupied'] = second_occupied
        second_complete['system'] = dir_name + str(int(i))
        second_complete['interface'] = '2nd'
        second_complete['ligand_type'] = ligand_type
        second_complete['ligand_state'] = ligand_state

        all_complete = pd.concat([all_complete, second_complete])

all_complete.to_csv('{}_{}_lipids.csv'.format(ligand_type, ligand_state))