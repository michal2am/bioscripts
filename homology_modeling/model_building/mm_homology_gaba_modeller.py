# python 2
# custom class for modeller GABA models
# michaladammichalowski@gmail.com
# ? - creation
# 14.12.15 - refactor

from modeller.automodel import *
from modeller import *


class GABAModel(automodel):
    def __init__(self, env, alnfile, knowns, sequence, assess_methods, segments, start_res):
        automodel.__init__(self, alnfile=alnfile, knowns=knowns, sequence=sequence,
                           assess_methods=assess_methods, env=env,)
        self.segments_names = segments
        self.start_res = start_res

    def special_patches(self, aln):
        segment_names = ['A', 'B', 'C', 'D', 'E']
        start_res = [8, 10, 8, 10, 23]
        self.rename_segments(segment_ids=segment_names, renumber_residues=start_res)
        self.patch(residue_type='DISU', residues=(self.residues['136:A'], self.residues['150:A']))
        self.patch(residue_type='DISU', residues=(self.residues['136:C'], self.residues['150:C']))
        self.patch(residue_type='DISU', residues=(self.residues['138:B'], self.residues['152:B']))
        self.patch(residue_type='DISU', residues=(self.residues['138:D'], self.residues['152:D']))
        self.patch(residue_type='DISU', residues=(self.residues['233:B'], self.residues['292:B']))
        self.patch(residue_type='DISU', residues=(self.residues['233:D'], self.residues['298:D']))
        self.patch(residue_type='DISU', residues=(self.residues['151:E'], self.residues['165:E']))
        self.patch(residue_type='DISU', residues=(self.residues['244:E'], self.residues['303:E']))
