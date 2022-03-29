import MDAnalysis
import numpy as np
import pandas as pd
from tqdm import tqdm
from MDAnalysis.analysis import distances

pdb = 'step5_input_chains.pdb'
xtc = 'all_runsteps.xtc'

u = MDAnalysis.Universe(pdb, xtc)

selection = 'segid {} and resid {} and name {}'
phe289 = selection.format('A', 289, 'CA')
met236 = selection.format('B', 236, 'CA')
asn265 = selection.format('A', 265, 'CA')

phe289 = u.select_atoms(phe289)
met236 = u.select_atoms(met236)
asn265 = u.select_atoms(asn265)

distances_m3m1 = []
distances_m3m2 = []

for ts in tqdm(u.trajectory):
    dist_arr1 = distances.distance_array(phe289.positions,  # reference
                                        met236.positions,  # configuration
                                        box=u.dimensions)
    distances_m3m1.append(dist_arr1[0][0])
    dist_arr2 = distances.distance_array(phe289.positions,  # reference
                                        asn265.positions,  # configuration
                                        box=u.dimensions)
    distances_m3m2.append(dist_arr2[0][0])

dists = pd.DataFrame()
dists['M3M1'] = distances_m3m1
dists['M3M2'] = distances_m3m2
dists.to_csv('dists.csv')
