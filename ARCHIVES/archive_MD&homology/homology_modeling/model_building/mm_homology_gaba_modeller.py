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
        self.patch(residue_type='DISU', residues=(self.residues['136:B'], self.residues['150:B']))
        self.patch(residue_type='DISU', residues=(self.residues['136:E'], self.residues['150:E']))
        self.patch(residue_type='DISU', residues=(self.residues['139:A'], self.residues['153:A']))
        self.patch(residue_type='DISU', residues=(self.residues['139:D'], self.residues['153:D']))
        # self.patch(residue_type='DISU', residues=(self.residues['233:B'], self.residues['292:B']))
        # self.patch(residue_type='DISU', residues=(self.residues['233:D'], self.residues['298:D']))
        self.patch(residue_type='DISU', residues=(self.residues['151:C'], self.residues['165:C']))
        # self.patch(residue_type='DISU', residues=(self.residues['244:E'], self.residues['303:E']))

        # self.patch(residue_type='DISU', residues=(self.residues['31:E'], self.residues['15:D']))
        # self.patch(residue_type='DISU', residues=(self.residues['31:B'], self.residues['15:A']))

        # remove bottom


    def special_restraints(self, aln):
        rsr = self.restraints
        at=self.atoms
        rsr.add(secondary_structure.alpha(self.residue_range('306:A', '311:A'))) #ashift 305:310 bshift 306:311 ab 311
        # rsr.add(secondary_structure.alpha(self.residue_range('315:A', '347:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('314:A', '322:A')))
        rsr.add(secondary_structure.alpha(self.residue_range('306:D', '311:D')))
        # rsr.add(secondary_structure.alpha(self.residue_range('315:D', '347:D')))
        rsr.add(secondary_structure.alpha(self.residue_range('314:D', '322:D')))

        rsr.add(secondary_structure.alpha(self.residue_range('300:B', '305:B'))) #ashift 300:305 bshift 301:306 ab 305
        # rsr.add(secondary_structure.alpha(self.residue_range('310:B', '340:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('310:B', '317:B')))
        rsr.add(secondary_structure.alpha(self.residue_range('300:E', '305:E')))
        # rsr.add(secondary_structure.alpha(self.residue_range('310:E', '340:E')))
        rsr.add(secondary_structure.alpha(self.residue_range('310:E', '317:E')))

        rsr.add(secondary_structure.alpha(self.residue_range('315:C', '320:C'))) #ashift 315:320 bshift 316:321 ab 320
        # rsr.add(secondary_structure.alpha(self.residue_range('325:C', '357:C')))
        rsr.add(secondary_structure.alpha(self.residue_range('325:C', '332:C')))


        '''
        chain_id = ['A', 'A', 'A', 'A', 'D', 'D', 'D', 'D',
                    'B', 'B', 'B', 'B', 'E', 'E', 'E', 'E',
                    'C', 'C', 'C', 'C',
                    'A', 'A', 'A',
                    'D', 'D', 'D']
        res_id = [311, 312, 313, 314, 311, 312, 313, 314,
                  306, 307, 308, 309, 306, 307, 308, 309,
                  321, 322, 323, 324,
                  173, 174, 175,
                  173, 174, 175]
        phi = [-1.23, -2.05, -0.56, -2.45, -1.23, -2.05, -0.56, -2.45,
               -0.99, -1.21, -1.36, -2.31, -0.99, -1.21, -1.36, -2.31,
               -1.23, -2.05, -0.56, -2.4,
               0.81, 1.06, 0.73,
               1.08, 0.98, -0.66]
        psi = [-0.27, -0.64, -1.03, 1.20, -0.27, -0.64, -1.03, 1.20,
               -0.39, -0.96, -0.54, 1.16, -0.39, -0.96, -0.54, 1.16,
               -0.27, -0.64, -1.03, 1.20,
               1.15, 1.44, -0.42,
               1.08, 1.47, -0.50]
        '''

        '''
        chain_id = ['A', 'A', 'A', 'A', 'D', 'D', 'D', 'D',
                    'B', 'B', 'B', 'B', 'E', 'E', 'E', 'E',
                    'C', 'C', 'C', 'C',
                    'A', 'A',
                    'D', 'D']
        res_id = [311, 312, 313, 314, 311, 312, 313, 314,
                  306, 307, 308, 309, 306, 307, 308, 309,
                  321, 322, 323, 324,
                  173, 174,
                  173, 174]
        phi = [-1.23, -2.05, -0.56, -2.45, -1.23, -2.05, -0.56, -2.45,
               -0.99, -1.21, -1.36, -2.31, -0.99, -1.21, -1.36, -2.31,
               -1.23, -2.05, -0.56, -2.4,
               -2.41, -1.02,
               -2.02, -1.04]
        psi = [-0.27, -0.64, -1.03, 1.20, -0.27, -0.64, -1.03, 1.20,
               -0.39, -0.96, -0.54, 1.16, -0.39, -0.96, -0.54, 1.16,
               -0.27, -0.64, -1.03, 1.20,
               -2.86, 2.44,
               -2.91, 2.49]
        '''

        chain_id = ['A', 'A', 'A', 'A', 'D', 'D', 'D', 'D',
                    'B', 'B', 'B', 'B', 'E', 'E', 'E', 'E',
                    'C', 'C', 'C', 'C',
                    'A',
                    'D']
        res_id = [311, 312, 313, 314, 311, 312, 313, 314,
                  306, 307, 308, 309, 306, 307, 308, 309,
                  321, 322, 323, 324,
                  174,
                  174]
        phi = [-1.23, -2.05, -0.56, -2.45, -1.23, -2.05, -0.56, -2.45,
               -0.99, -1.21, -1.36, -2.31, -0.99, -1.21, -1.36, -2.31,
               -1.23, -2.05, -0.56, -2.4,
               -1.02,
               -1.04]
        psi = [-0.27, -0.64, -1.03, 1.20, -0.27, -0.64, -1.03, 1.20,
               -0.39, -0.96, -0.54, 1.16, -0.39, -0.96, -0.54, 1.16,
               -0.27, -0.64, -1.03, 1.20,
               2.44,
               2.49]

        for chain_id, res_id, phi, psi in zip(chain_id, res_id, phi, psi):

            rsr.add(forms.gaussian(group=physical.phi_dihedral,
                                   feature=features.dihedral([at['C:' + str(res_id - 1) + ':' + chain_id],
                                                              at['N:' + str(res_id) + ':' + chain_id],
                                                              at['CA:' + str(res_id) + ':' + chain_id],
                                                              at['C:' + str(res_id) + ':' + chain_id]]),
                                   mean=phi, stdev=0.1))
            rsr.add(forms.gaussian(group=physical.psi_dihedral,
                                   feature=features.dihedral([at['N:' + str(res_id) + ':' + chain_id],
                                                              at['CA:' + str(res_id) + ':' + chain_id],
                                                              at['C:' + str(res_id) + ':' + chain_id],
                                                              at['N:' + str(res_id + 1) + ':' + chain_id]]),
                                   mean=psi, stdev=0.1))

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
