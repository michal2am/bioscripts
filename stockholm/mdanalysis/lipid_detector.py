import numpy as np
import MDAnalysis
import MDAnalysis.analysis.align
from tqdm import tqdm
import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--pdb") # 'step5_input_chains.pdb'
parser.add_argument("--xtc") # 'all_runsteps.xtc'
parser.add_argument("--dir") # 'sys'
parser.add_argument("--l_type") # 'etomidate'
parser.add_argument("--l_state") # 'apo'

args = parser.parse_args()

pdb = args.pdb
xtc_f = args.xtc
dir_name = args.dir
ligand_type = args.l_type
ligand_state = args.l_state


class SubunitInterface:

    def __init__(self, interface_name, atom_selection, system, l_type, l_state):
        self.interface_name = interface_name
        self.atom_selection = atom_selection
        self.system = system
        self.ligand_type = l_type
        self.lignad_state = l_state

        self.lipids = U.select_atoms('((resname POPC) and sphzone 5.0 {})'.format(self.atom_selection), updating=True)
        self.lipids_contacts = []

    def lipid_data(self):
        lipids_data = pd.DataFrame()
        lipids_data['occupied'] = self.lipids_contacts
        lipids_data['system'] = self.system
        lipids_data['interface'] = self.interface_name
        lipids_data['ligand_type'] = self.ligand_type
        lipids_data['ligand_state'] = self.lignad_state
        return lipids_data


interface_configs = [['1st_ba', '(segid A and resid 289) or (segid A and resid 286) or (segid A and resid 265) or (segid C and resid 232) or (segid C and resid 233)'],
                     ['2nd_ba', '(segid B and resid 289) or (segid B and resid 286) or (segid B and resid 265) or (segid D and resid 232) or (segid D and resid 233)']
                    ]
interfaces = []

for i in range(1, 2):

    print('Going to sys{}'.format(i))
    # xtc = dir_name + str(int(i)) + '/gromacs/' + xtc_f
    xtc = dir_name + str(int(i)) + '/' + xtc_f

    if os.path.isfile(xtc):

        U = MDAnalysis.Universe(pdb, xtc)

        for interface_config in interface_configs:
            interfaces.append(SubunitInterface(interface_config[0], interface_config[1],
                                               'MD' + str(i), ligand_type, ligand_state))

        for ts in tqdm(U.trajectory):
            for interface in interfaces:
                interface.lipids_contacts.append(0) if len(interface.lipids) == 0 else interface.lipids_contacts.append(1)

all_complete = pd.DataFrame()

for interface in interfaces:
    single_complete = interface.lipid_data()
    all_complete = pd.concat([all_complete, single_complete])

all_complete.index.name = 'frame'
all_complete.to_csv('{}_{}_lipids.csv'.format(ligand_type, ligand_state))