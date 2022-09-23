# TODO: needs small refactor like gate_resid_dihe

import MDAnalysis
from tqdm import tqdm
import pandas as pd
import argparse
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument("--pdb") # full path to pdb
parser.add_argument("--xtc") # xtc file name
parser.add_argument("--dir") # full path to the dir with dir prefix /xxx/yyy/MD
parser.add_argument("--num", type=int) # replica number
parser.add_argument("--l_type") # ligand name
parser.add_argument("--l_state") # apo/holo

args = parser.parse_args()

pdb = args.pdb
xtc_f = args.xtc
dir_name = args.dir
dir_num = args.num
xtc = dir_name + str(int(dir_num)) + '/' + xtc_f
# xtc = dir_name + str(int(i)) + '/gromacs/' + xtc_f
ligand_type = args.l_type
ligand_state = args.l_state


class SubunitInterface:
    """
    object to store the data, nothing using the mda universe as it is
    """

    def __init__(self, interface_name,  system, l_type, l_state, gate_res_selection):
        self.interface_name = interface_name
        self.system = system
        self.ligand_type = l_type
        self.lignad_state = l_state

        self.prim_res = gate_res_selection[0]
        self.comp_res = gate_res_selection[1]
        self.gate_dist = []

    def gate_resi_data(self):
        data = pd.DataFrame()
        data['distance'] = self.gate_dist
        data['system'] = self.system
        data['interface'] = self.interface_name
        data['ligand_type'] = self.ligand_type
        data['ligand_state'] = self.lignad_state
        return data

order = 'ACBDE'

if order == 'ABCDE':
    # just for bicuculline A-B-C-D-E notation
    interface_configs = [['1st_beta/alpha', ['segid A and name SD and resid 286 and protein', 'segid B and name CG and resid 232 and protein']],
                         ['2nd_beta/alpha', ['segid C and name SD and resid 286 and protein', 'segid D and name CG and resid 232 and protein']],
                         ]
if order == 'ACBDE':
    # just for others A-C-B-D-E notation
    interface_configs = [['1st_beta/alpha', ['segid A and name SD and resid 286 and protein', 'segid C and name CG and resid 232 and protein']],
                         ['2nd_beta/alpha', ['segid B and name SD and resid 286 and protein', 'segid D and name CG and resid 232 and protein']],
                         ]

# mda universe
universe = MDAnalysis.Universe(pdb, xtc)
# container for obtained data
all_interfaces = []

# create data container for each subunit interface
for interface_config in interface_configs:
    all_interfaces.append(SubunitInterface(interface_config[0], 'MD' + str(dir_num), ligand_type, ligand_state, interface_config[1],))

# fill each container with data while going through trajectory
for ts in tqdm(universe.trajectory):
    for interface in all_interfaces:
        prim = universe.select_atoms(interface.prim_res)
        comp = universe.select_atoms(interface.comp_res)
        interface.gate_dist.append(np.linalg.norm(prim.center_of_mass() - comp.center_of_mass()))


# gather and concatenate data from each subunit interface
all_complete = pd.DataFrame()

for interface in all_interfaces:
    all_complete = pd.concat([all_complete, interface.gate_resi_data()])

all_complete.index.name = 'frame'
all_complete.to_csv('{}_{}_MD{}_gate.csv'.format(ligand_type, ligand_state, dir_num))