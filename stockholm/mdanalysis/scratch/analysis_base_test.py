import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
import numpy as np


def lipid_in(lipid_list):
    occupied = 0 if len(lipid_list) == 0 else 1
    return occupied

class LipidContacts(AnalysisBase):

    def __init__(self, atomgroup, verbose=True):
        trajectory = atomgroup.universe.trajectory
        super().__init__(trajectory, verbose=True)

        self.atomgroup = atomgroup

    def _prepare(self):
        self.results.contacts = np.zeros((self.n_frames, 1))

    def _single_frame(self):
        contacts = lipid_in(self.atomgroup)
        self.results.contacts[self._frame_index] = contacts

    def _conclude(self):
        self.contacts_sum = np.sum(self.results.contacts)/len(self.results.contacts)

pdb = '/mnt/cephfs/projects/2022033000_GABA_mm_internship/as_retraj/Prop_Apo/view.pdb'
xtc = '/mnt/cephfs/projects/2022033000_GABA_mm_internship/as_retraj/Prop_Apo/MD1/view.xtc'


u = mda.Universe(pdb, xtc)

lipids_ba1st = u.select_atoms('((resname POPC) and sphzone 5.0 {})'.format
                              ('((segid A and resid 289) or (segid A and resid 286) or (segid A and resid 265) '
                               'or (segid C and resid 232))'), updating=True)

test_run = LipidContacts(lipids_ba1st).run()
print(test_run.contacts_sum)
