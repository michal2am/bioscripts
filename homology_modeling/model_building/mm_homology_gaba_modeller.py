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
        self.rename_segments(segment_ids=self.segments_names, renumber_residues=self.start_res)
        self.patch(residue_type='DISU', residues=(self.residues['136:A'], self.residues['150:A']))
        self.patch(residue_type='DISU', residues=(self.residues['136:C'], self.residues['150:C']))
        self.patch(residue_type='DISU', residues=(self.residues['139:B'], self.residues['153:B']))
        self.patch(residue_type='DISU', residues=(self.residues['139:D'], self.residues['153:D']))
        # self.patch(residue_type='DISU', residues=(self.residues['233:B'], self.residues['292:B']))
        # self.patch(residue_type='DISU', residues=(self.residues['233:D'], self.residues['298:D']))
        self.patch(residue_type='DISU', residues=(self.residues['151:E'], self.residues['165:E']))
        # self.patch(residue_type='DISU', residues=(self.residues['244:E'], self.residues['303:E']))

'''
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
        self.rename_segments(segment_ids=self.segments_names, renumber_residues=self.start_res)

    #def select_atoms(self):
    #    sel = selection(self.residue_range('193:A', '209:A'))
    #    sel.add(self.residue_range('193:C', '209:C'))
    #    sel.add(self.residue_range('197:B', '213:B'))
    #    sel.add(self.residue_range('197:D', '213:D'))
    #    sel.add(self.residue_range('208:E', '224:E'))
    #    return sel

    def special_restraints(self, aln):
        rsr = self.restraints

        rsr.add(secondary_structure.strand(self.residue_range('193:A', '200:A')))
        rsr.add(secondary_structure.strand(self.residue_range('203:A', '209:A')))
        rsr.add(secondary_structure.strand(self.residue_range('193:C', '200:C')))
        rsr.add(secondary_structure.strand(self.residue_range('203:C', '209:C')))

        rsr.add(secondary_structure.strand(self.residue_range('197:B', '204:B')))
        rsr.add(secondary_structure.strand(self.residue_range('207:B', '213:B')))
        rsr.add(secondary_structure.strand(self.residue_range('197:D', '204:D')))
        rsr.add(secondary_structure.strand(self.residue_range('207:D', '213:D')))

        rsr.add(secondary_structure.strand(self.residue_range('208:E', '215:E')))
        rsr.add(secondary_structure.strand(self.residue_range('218:E', '224:E')))
'''