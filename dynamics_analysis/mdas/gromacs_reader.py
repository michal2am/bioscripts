import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align, rms

'''
mobile = mda.Universe('step5_input.gro', 'step7_1.trr')
ref = mda.Universe('step5_input.gro', 'step7_1.trr')

mobile.trajectory[-1]
ref.trajectory[0]

mobile_ca=mobile.select_atoms('name CA')
ref_ca = ref.select_atoms('name CA')
print(rms.rmsd(mobile_ca.positions, ref_ca.positions, superposition=False))

aligner = align.AlignTraj(mobile, ref, select='name CA', in_memory=True).run()

mobile.trajectory[-1]
ref.trajectory[0]

mobile_ca=mobile.select_atoms('name CA')
ref_ca = ref.select_atoms('name CA')
print(rms.rmsd(mobile_ca.positions, ref_ca.positions, superposition=False))

to_save = mobile.select_atoms('all')
to_save.write('test.trr', frames='all')
to_save.write('test.xtc', frames='all')
'''

u = mda.Universe('step7_1.tpr', 'all_runsteps.xtc', in_memory=False)
protein = u.select_atoms('protein')
for ts in u.trajectory:
    print(ts)
    protein.unwrap(compound='fragments')
for ts in u.trajectory:
    protein_center = protein.center_of_mass(pbc=True)
    dim = ts.triclinic_dimensions
    box_center = np.sum(dim, axis=0) / 2
    u.atoms.translate(box_center - protein_center)
not_protein = u.select_atoms('not protein')

for ts in u.trajectory:
    not_protein.wrap(compound='residues')

to_save = u.select_atoms('all')
to_save.write('test.xtc', frames='all')