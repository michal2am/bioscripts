from modeller import *
from modeller.automodel import *   

class GabaModel(automodel):
    def special_patches(self, aln):
        self.rename_segments(segment_ids=['D', 'C', 'B', 'A', 'E'], renumber_residues=[12, 10, 12, 10, 25])
