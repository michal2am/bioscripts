from modeller import *
from modeller.automodel import *


class GABArefine(loopmodel):

    def __init__(self, env, inimodel, sequence, location):
        loopmodel.__init__(self, env, inimodel=inimodel, sequence=sequence)
        self.location = location

    def select_loop_atoms(self):

        if self.location == "b4b5inner":

            sel = selection(self.residue_range('100:A', '106:A'))
            sel.add(self.residue_range('100:C', '106:C'))
            sel.add(self.residue_range('102:B', '108:B'))
            sel.add(self.residue_range('102:D', '108:D'))
            sel.add(self.residue_range('115:E', '121:E'))

        if self.location == "preLC":

            sel = selection(self.residue_range('195:A', '203:A'))
            sel.add(self.residue_range('194:C', '201:C'))
            sel.add(self.residue_range('199:B', '205:B'))
            sel.add(self.residue_range('199:D', '205:D'))
            sel.add(self.residue_range('210:E', '217:E'))

        return sel




