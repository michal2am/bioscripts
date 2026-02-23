import MDAnalysis
import numpy as np
from tqdm import tqdm


pdb = 'step5_input_chains.pdb'
xtc = 'all_runsteps.xtc'

U = MDAnalysis.Universe(pdb,xtc)
dists = np.zeros((U.trajectory.n_frames,2), dtype = float)
top_M2_sel = ['segid A and name CA and resid 266:271', 'segid C and name CA and resid 266:271']
tmd = U.select_atoms('name CA and ((segid A and resid 218:400) or (segid B and resid 223:400) or (segid C and resid 218:400) or (segid D and resid 223:400) or (segid E and resid 233:400))')

top_M2 = []
for j in top_M2_sel:
    top_M2.append(U.select_atoms(j))

for ts in tqdm(U.trajectory):
    for j in range(len(top_M2)):
        dists[ts.frame,j] = np.linalg.norm(tmd.center_of_mass() - top_M2[j].center_of_mass())

np.savetxt('M2_spread.dat', dists, newline='\n', delimiter='\t')

