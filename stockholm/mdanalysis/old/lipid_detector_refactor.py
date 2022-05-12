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
parser.add_argument("--dir") # '/xxxx/yyyy/sys'
parser.add_argument("--num", type=int) # '4'
parser.add_argument("--l_type") # 'etomidate'
parser.add_argument("--l_state") # 'apo'

args = parser.parse_args()

pdb = args.pdb
xtc_f = args.xtc
dir_name = args.dir
dir_num = args.num
ligand_type = args.l_type
ligand_state = args.l_state


class CompleteSystem:

    def __init__(self, pdb, traj, system, l_type, l_state):
        self.universe = MDAnalysis.Universe(pdb, traj)
        self.system = system
        self.ligand_type = l_type
        self.lignad_state = l_state


class SubunitInterfaceLipids:

    def __init__(self, complete_system, interface_name, lipid_contact_selection):
        self.interface_name = interface_name
        self.complete_system = complete_system

        self.lipid_contact_selection = lipid_contact_selection
        self.lipids = self.complete_system.universe.select_atoms('((resname POPC) and sphzone 5.0 {})'.format(self.lipid_contact_selection), updating=True)
        self.lipids_contacts = []

    def lipid_data(self):
        lipids_data = pd.DataFrame()
        lipids_data['occupied'] = self.lipids_contacts
        lipids_data['system'] = self.system
        lipids_data['interface'] = self.interface_name
        lipids_data['ligand_type'] = self.ligand_type
        lipids_data['ligand_state'] = self.lignad_state
        return lipids_data

# proline included:
# '1st_ba', '(segid A and resid 289) or (segid A and resid 286) or (segid A and resid 265) or (segid C and resid 232) or (segid C and resid 233)'
# '2nd_ba', '(segid B and resid 289) or (segid B and resid 286) or (segid B and resid 265) or (segid D and resid 232) or (segid D and resid 233)'

# [interface name, atom selection for lipid contacts, XXX]
interface_configs = [['1st_beta/alpha', '(segid A and resid 289) or (segid A and resid 286) or (segid A and resid 265) or (segid C and resid 232)'],
                     ['2nd_beta/alpha', '(segid B and resid 289) or (segid B and resid 286) or (segid B and resid 265) or (segid D and resid 232)'],
                     ['alpha/beta', '(segid C and resid 294) or (segid C and resid 291) or (segid C and resid 270) or (segid B and resid 227)'],
                     ['alpha/gamma', '(segid D and resid 294) or (segid D and resid 291) or (segid D and resid 270) or (segid E and resid 242)'],
                     ['gamma/beta', '(segid E and resid 304) or (segid E and resid 301) or (segid E and resid 280) or (segid A and resid 227)']]


def check_interfaces(i):

    current_sys = 'sys{}'.format(i)
    xtc = dir_name + str(int(i)) + '/' + xtc_f

    complete_system = CompleteSystem(pdb, xtc, 'MD' + str(i), ligand_type, ligand_state)
    interfaces_system = []

    for interface_config in interface_configs:
        interfaces_system.append(SubunitInterfaceLipids(complete_system, interface_config[0], interface_config[1]))

        for ts in tqdm(U.trajectory):
            for interface in interfaces_all_systems[current_sys]:
                interface.lipids_distances.append(0) if len(interface.lipids) == 0 else interface.lipids_distances.append(1)

for i in range(1, dir_num + 1):
    check_interfaces(i)

all_complete = pd.DataFrame()

for system, interfaces in interfaces_all_systems.items():
    for interface in interfaces:
        single_complete = interface.lipid_data()
        all_complete = pd.concat([all_complete, single_complete])

all_complete.index.name = 'frame'
all_complete.to_csv('{}_{}_lipids.csv'.format(ligand_type, ligand_state))