import MDAnalysis
from tqdm import tqdm
import pandas as pd
import argparse

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

    def __init__(self, interface_name,  system, l_type, l_state, lipid_contact_selection=None):
        self.interface_name = interface_name
        self.system = system
        self.ligand_type = l_type
        self.lignad_state = l_state

        # only for lipids:
        self.lipids_selection = '((resname POPC) and sphzone 5.0 {})'.format(lipid_contact_selection)
        self.lipids_contacts = []

    def lipid_data(self):
        lipids_data = pd.DataFrame()
        lipids_data['occupied'] = self.lipids_contacts
        lipids_data['system'] = self.system
        lipids_data['interface'] = self.interface_name
        lipids_data['ligand_type'] = self.ligand_type
        lipids_data['ligand_state'] = self.lignad_state
        return lipids_data

'''
interface_configs = [['1st_beta/alpha', '((segid A and resid 289) or (segid A and resid 286) or (segid A and resid 265) or (segid C and resid 232))'],
                     ['2nd_beta/alpha', '((segid B and resid 289) or (segid B and resid 286) or (segid B and resid 265) or (segid D and resid 232))'],
                     ['alpha/beta', '((segid C and resid 294) or (segid C and resid 291) or (segid C and resid 270) or (segid B and resid 227))'],
                     ['alpha/gamma', '((segid D and resid 294) or (segid D and resid 291) or (segid D and resid 270) or (segid E and resid 242))'],
                     ['gamma/beta', '((segid E and resid 304) or (segid E and resid 301) or (segid E and resid 280) or (segid A and resid 227))']]
'''
interface_configs = [['1st_beta/alpha', '((segid A and resid 289) or (segid A and resid 286) or (segid A and resid 265) or (segid B and resid 232))'],
                     ['2nd_beta/alpha', '((segid C and resid 289) or (segid C and resid 286) or (segid C and resid 265) or (segid D and resid 232))'],
                     ['alpha/beta', '((segid B and resid 294) or (segid B and resid 291) or (segid B and resid 270) or (segid C and resid 227))'],
                     ['alpha/gamma', '((segid D and resid 294) or (segid D and resid 291) or (segid D and resid 270) or (segid E and resid 242))'],
                     ['gamma/beta', '((segid E and resid 304) or (segid E and resid 301) or (segid E and resid 280) or (segid A and resid 227))']]

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
        lipid_atoms = universe.select_atoms(interface.lipids_selection)
        interface.lipids_contacts.append(0) if len(lipid_atoms) == 0 else interface.lipids_contacts.append(1)

# gather and concatenate data from each subunit interface
all_complete = pd.DataFrame()

for interface in all_interfaces:
    all_complete = pd.concat([all_complete, interface.lipid_data()])

all_complete.index.name = 'frame'
all_complete.to_csv('{}_{}_MD{}_lipids.csv'.format(ligand_type, ligand_state, dir_num))
