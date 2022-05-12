import MDAnalysis
from tqdm import tqdm
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
#python ~/repos/bioscripts/stockholm/mdanalysis/lipid_detector.py --pdb /mnt/my_cephfs/as_retraj/Phenobarb/sim/phenobarb_holo_protligpopcion.pdb --dir /mnt/my_cephfs/as_retraj/Phenobarb/sim/MD --l_type phenobarbital --l_state holo --ord ABCDE --xtc phenobarb_holo_protligpopcion_MD2.xtc --num 2 --rad 6

parser.add_argument("--pdb")    # full path to pdb
parser.add_argument("--ord")    # chain id order
parser.add_argument("--xtc")    # xtc file name
parser.add_argument("--dir")    # full path to the dir with dir prefix /xxx/yyy/MD
parser.add_argument("--num", type=int)  # replica number
parser.add_argument("--l_type")     # ligand name
parser.add_argument("--l_state")    # apo/holo
parser.add_argument("--rad")    # sphere detector radius

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
radius = args.rad


class SubunitInterface:
    """
    object to store the data, nothing using the mda universe as it is
    """

    def __init__(self, interface_name,  system, l_type, l_state, rad, lipid_contact_selection=None):
        self.interface_name = interface_name
        self.system = system
        self.ligand_type = l_type
        self.lignad_state = l_state

        # only for lipids:
        self.det_rad = rad
        self.lipids_selection = '((resname POPC) and sphzone {} {})'.format(self.det_rad, lipid_contact_selection)
        self.lipids_contacts = []

    def data(self):
        data = pd.DataFrame()
        data['occupied'] = self.lipids_contacts
        data['system'] = self.system
        data['interface'] = self.interface_name
        data['ligand_type'] = self.ligand_type
        data['ligand_state'] = self.lignad_state
        data['detection_radius'] = self.det_rad
        return data


if chain_order == 'ACBDE':
    interface_configs = [['1st_beta/alpha', '(name CA and ((segid A and resid 289) or (segid A and resid 265) or (segid C and resid 233)))'],
                         ['2nd_beta/alpha', '(name CA and ((segid B and resid 289) or (segid B and resid 265) or (segid D and resid 233)))'],

                         ['alpha/beta', '(name CA and ((segid C and resid 294) or (segid C and resid 270) or (segid B and resid 228)))'],
                         ['alpha/gamma', '(name CA and ((segid D and resid 294) or (segid D and resid 270) or (segid E and resid 243)))'],

                         ['gamma/beta', '(name CA and ((segid E and resid 304) or (segid E and resid 280) or (segid A and resid 228)))']]
    '''
    # original:
    interface_configs = [['1st_beta/alpha', '((segid A and resid 289) or (segid A and resid 286) or '
                                            '(segid A and resid 265) or (segid C and resid 232))'],
                         ['2nd_beta/alpha', '((segid B and resid 289) or (segid B and resid 286) or '
                                            '(segid B and resid 265) or (segid D and resid 232))'],
                         ['alpha/beta', '((segid C and resid 294) or (segid C and resid 291) or '
                                        '(segid C and resid 270) or (segid B and resid 227))'],
                         ['alpha/gamma', '((segid D and resid 294) or (segid D and resid 291) or '
                                         '(segid D and resid 270) or (segid E and resid 242))'],
                         ['gamma/beta', '((segid E and resid 304) or (segid E and resid 301) or '
                                        '(segid E and resid 280) or (segid A and resid 227))']]
    '''

if chain_order == 'ABCDE':
    interface_configs = [['1st_beta/alpha', '(name CA and ((segid A and resid 289) or (segid A and resid 265) or (segid B and resid 233)))'],
                         ['2nd_beta/alpha', '(name CA and ((segid C and resid 289) or (segid C and resid 265) or (segid D and resid 233)))'],

                         ['alpha/beta', '(name CA and ((segid B and resid 294) or (segid B and resid 270) or (segid C and resid 228)))'],
                         ['alpha/gamma', '(name CA and ((segid D and resid 294) or (segid D and resid 270) or (segid E and resid 243)))'],

                         ['gamma/beta', '(name CA and ((segid E and resid 304) or (segid E and resid 280) or (segid A and resid 228)))']]
    '''
    # original:
    interface_configs = [['1st_beta/alpha', '((segid A and resid 289) or (segid A and resid 286) or '
                                            '(segid A and resid 265) or (segid B and resid 232))'],
                         ['2nd_beta/alpha', '((segid C and resid 289) or (segid C and resid 286) or '
                                            '(segid C and resid 265) or (segid D and resid 232))'],
                         ['alpha/beta', '((segid B and resid 294) or (segid B and resid 291) or '
                                        '(segid B and resid 270) or (segid C and resid 227))'],
                         ['alpha/gamma', '((segid D and resid 294) or (segid D and resid 291) or '
                                         '(segid D and resid 270) or (segid E and resid 242))'],
                         ['gamma/beta', '((segid E and resid 304) or (segid E and resid 301) or '
                                        '(segid E and resid 280) or (segid A and resid 227))']]
    '''

# mda universe
universe = MDAnalysis.Universe(pdb, xtc)

# create data container for each subunit interface
all_interfaces = [SubunitInterface(interface_config[0], 'MD' + str(dir_num), ligand_type, ligand_state,
                                   radius, interface_config[1])
                  for interface_config in interface_configs]

# fill each container with data while going through trajectory
for ts in tqdm(universe.trajectory):
    for interface in all_interfaces:
        lipid_atoms = universe.select_atoms(interface.lipids_selection)
        interface.lipids_contacts.append(0) if len(lipid_atoms) == 0 else interface.lipids_contacts.append(1)

# gather and concatenate data from each subunit interface
all_complete = pd.DataFrame()

for interface in all_interfaces:
    all_complete = pd.concat([all_complete, interface.data()])

all_complete.index.name = 'frame'
all_complete.to_csv('{}_{}_{}_MD{}_lipids.csv'.format(ligand_type, ligand_state, radius, dir_num))
