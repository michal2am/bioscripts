# same as phe_helices, but according to AS scripts for reproducibility and performance
# works on pre-processed xtc and corrected charmm-gui input .pdb

import numpy as np
import MDAnalysis
import MDAnalysis.analysis.align
from tqdm import tqdm
import os

pdb = 'step5_input_chains.pdb'
xtc_f = 'all_runsteps.xtc'

for i in range(1, 5):

    print('Going to sys{}'.format(i))

    xtc = 'sys' + str(int(i)) + '/gromacs/' + xtc_f
    print(xtc)
    if os.path.isfile(xtc):

        U = MDAnalysis.Universe(pdb, xtc)

        M3_1 = U.select_atoms('segid A and name CA and resid 286 and protein')
        M3_2 = U.select_atoms('segid B and name CA and resid 286 and protein')

        M1_1 = U.select_atoms('segid C and name CA and resid 232 and protein')
        M1_2 = U.select_atoms('segid D and name CA and resid 232 and protein')

        dist = np.zeros((U.trajectory.n_frames, 2), dtype=float)

        for ts in tqdm(U.trajectory):
            dist[ts.frame,0] = np.linalg.norm(M3_1.center_of_mass() - M1_1.center_of_mass())
            dist[ts.frame,1] = np.linalg.norm(M3_2.center_of_mass() - M1_2.center_of_mass())

        out_file = 'sys' + str(int(i)) + '/m3_m1.dat'
        np.savetxt(out_file, dist, newline = '\n', delimiter = '\t')



