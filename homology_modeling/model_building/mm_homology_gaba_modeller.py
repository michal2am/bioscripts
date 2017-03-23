# python 2
# custom class for modeller GABA models
# michaladammichalowski@gmail.com
# ? - creation
# 14.12.15 - refactor

from modeller.automodel import *
from modeller import *


class GABAModel(automodel):
    """
    basic model with disulph bridges
    """
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


class GABAModelSelection(automodel):
    """
    template based refinement of LC region
    """
    def __init__(self, env, alnfile, knowns, sequence, assess_methods, segments, start_res):
        automodel.__init__(self, alnfile=alnfile, knowns=knowns, sequence=sequence,
                           assess_methods=assess_methods, env=env,)
        self.segments_names = segments
        self.start_res = start_res

    def special_patches(self, aln):
        segment_names = ['A', 'B', 'C', 'D', 'E']
        start_res = [8, 10, 8, 10, 23]
        self.rename_segments(segment_ids=segment_names, renumber_residues=start_res)

    def select_atoms(self):
        sel = selection(self.residue_range('195:A', '208:A'))
        sel.add(self.residue_range('195:C', '208:C'))
        sel.add(self.residue_range('199:B', '212:B'))
        sel.add(self.residue_range('199:D', '212:D'))
        sel.add(self.residue_range('210:E', '223:E'))
        return sel

    def special_restraints(self, aln):
        rsr = self.restraints

        rsr.add(secondary_structure.strand(self.residue_range('195:A', '199:A')))
        rsr.add(secondary_structure.strand(self.residue_range('204:A', '208:A')))
        rsr.add(secondary_structure.strand(self.residue_range('195:C', '199:C')))
        rsr.add(secondary_structure.strand(self.residue_range('204:C', '208:C')))

        rsr.add(secondary_structure.strand(self.residue_range('199:B', '203:B')))
        rsr.add(secondary_structure.strand(self.residue_range('208:B', '212:B')))
        rsr.add(secondary_structure.strand(self.residue_range('199:D', '203:D')))
        rsr.add(secondary_structure.strand(self.residue_range('208:D', '212:D')))

        rsr.add(secondary_structure.strand(self.residue_range('210:E', '214:E')))
        rsr.add(secondary_structure.strand(self.residue_range('219:E', '223:E')))
