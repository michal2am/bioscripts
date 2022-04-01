import MDAnalysis
import numpy as np
import pandas as pd
from tqdm import tqdm

# gmx trjconv -f step5_input.pdb -s step7_1.tpr -o step5_input_chains.pdb

pdb = 'step5_input_chains.pdb'
xtc = 'all_runsteps.xtc'

u = MDAnalysis.Universe(pdb, xtc)
psi_angles = []

for segid in ['A', 'C']:

    selection = 'segid {} and resname {} and resid {} and name {}'
    phe1 = selection.format(segid, 'PHE', 289, 'CG')
    phe2 = selection.format(segid, 'PHE', 289, 'CB')
    phe3 = selection.format(segid, 'PHE', 289, 'CA')
    phe4 = selection.format(segid, 'PHE', 289, 'C')

    psi = [phe1, phe2, phe3, phe4]
    psi_angle = sum([u.select_atoms(atom) for atom in psi])  # sum of Atoms creates an AtomGroup
    psi_angle = psi_angle.dihedral  # convert AtomGroup to Dihedral object
    psi_angles.append(psi_angle)

psi_list_a = []
psi_list_b = []

for ts in tqdm(u.trajectory):
    psi_list_a.append(psi_angles[0].value())
    psi_list_b.append(psi_angles[1].value())

psis = pd.DataFrame()
psis['segid_A'] = psi_list_a
psis['segid_C'] = psi_list_b
psis.to_csv('phe_psis.csv')