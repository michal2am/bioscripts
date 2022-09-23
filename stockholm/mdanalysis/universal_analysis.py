import MDAnalysis
from tqdm import tqdm
import pandas as pd
import argparse
import numpy as np

parser = argparse.ArgumentParser()
#python ~/repos/bioscripts/stockholm/mdanalysis/latch_resid_dihe.py --pdb /mnt/my_cephfs/as_retraj/Prop_Apo/propofol_apo_protligpopcion.pdb --dir /mnt/my_cephfs/as_retraj/Prop_Apo/MD  --l_type propofol --ord ACBDE --l_state apo --xtc propofol_apo_protligpopcion_MD1.xtc --num 1


parser.add_argument("--pdb")    # full path to pdb
parser.add_argument("--ord")    # chain id order
parser.add_argument("--xtc")    # xtc file name
parser.add_argument("--dir")    # full path to the dir with dir prefix /xxx/yyy/MD
parser.add_argument("--num", type=int)   # replica number
parser.add_argument("--l_type")     # ligand name
parser.add_argument("--l_state")    # apo/holo

args = parser.parse_args()

pdb = args.pdb
xtc_f = args.xtc
dir_name = args.dir
dir_num = args.num
xtc = dir_name + str(int(dir_num)) + '/' + xtc_f
# xtc = dir_name + str(int(i)) + '/gromacs/' + xtc_f
ligand_type = args.l_type
ligand_state = args.l_state
chain_order = args.ord

class SubunitInterface:
    """
    object to store the data, nothing using the mda universe as it is
    """

    def __init__(self, interface_name,  system, l_type, l_state, chains):
        self.interface_name = interface_name
        self.system = system
        self.ligand_type = l_type
        self.lignad_state = l_state
        self.princ_selection = chains[0]
        self.comp_selection = chains[1]

        self.data = {}

    def data(self):
        data = pd.DataFrame(self.data)
        data['system'] = self.system
        data['interface'] = self.interface_name
        data['ligand_type'] = self.ligand_type
        data['ligand_state'] = self.lignad_state
        return data


class Receptor:
    """

    """

    def __init__(self, interfaces):
        self.interfaces = interfaces

    def calculate_distance_per_interface(self, atoms_names, atoms_label):
        for interface in self.interfaces:
            interface.data[atoms_label] = []
        for ts in tqdm(universe.trajectory):
            for interface in self.interfaces:
                prim = universe.select_atoms(interface.princ_selection + atoms_names[0])
                comp = universe.select_atoms(interface.comp_selection + atoms_names[1])
                interface.data[atoms_label].append(np.linalg.norm(prim.center_of_mass() - comp.center_of_mass()))

if chain_order == 'ABCDE':
    interface_configs = [['1st_beta/alpha', ['segid A and ', 'segid B and ']],
                         ['2nd_beta/alpha', ['segid C and ', 'segid D and ']]]


# mda universe
universe = MDAnalysis.Universe(pdb, xtc)

# create data container for each subunit interface
receptor = Receptor([SubunitInterface(interface_config[0], 'MD' + str(dir_num), ligand_type, ligand_state, interface_config[1])
                     for interface_config in interface_configs])

receptor.calculate_distance_per_interface(['name CA and resid 200 and protein', 'name CA and resid 65 and protein'], 'loop C opening')
# fill each container with data while going through trajectory


# gather and concatenate data from each subunit interface

# TODO: those should use receptor object, not subunit
all_complete = pd.DataFrame()

for interface in receptor.interfaces:
    all_complete = pd.concat([all_complete, interface.data()])

all_complete.index.name = 'frame'
all_complete.to_csv('{}_{}_MD{}_distance.csv'.format(ligand_type, ligand_state, dir_num))