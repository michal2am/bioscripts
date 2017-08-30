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

        if self.location == 'postA1':

            sel = selection(self.residue_range('14:A', '20:A'))
            sel.add(self.residue_range('14:C', '20:C'))
            sel.add(self.residue_range('16:B', '22:B'))
            sel.add(self.residue_range('16:D', '22:D'))
            sel.add(self.residue_range('29:E', '35:E'))

        if self.location == 'preA2':

            sel = selection(self.residue_range('78:B', '84:B'))
            sel.add(self.residue_range('78:D', '84:D'))
            sel.add(self.residue_range('91:E', '97:E'))

        if self.location == 'LF':

            sel = selection(self.residue_range('175:B', '179:B'))
            #sel.add(self.residue_range('181:C', '187:C'))
            sel.add(self.residue_range('186:E', '190:E'))

        if self.location == 'c184':

            sel = selection(self.residue_range('182:C', '186:C'))

        if self.location == 'preB1B2':

            sel = selection(self.residue_range('49:B', '53:B'))
            sel.add(self.residue_range('49:D', '53:D'))
            sel.add(self.residue_range('62:E', '66:E'))

        if self.location == "preLC":

            sel = selection(self.residue_range('197:A', '199:A'))
            sel.add(self.residue_range('197:C', '199:C'))
            sel.add(self.residue_range('201:B', '203:B'))
            sel.add(self.residue_range('201:D', '203:D'))
            sel.add(self.residue_range('212:E', '214:E'))

        if self.location == "preLC2":

            sel = selection(self.residue_range('197:A', '199:A'))
            sel.add(self.residue_range('197:C', '199:C'))
            #sel.add(self.residue_range('201:B', '203:B'))
            #sel.add(self.residue_range('201:D', '203:D'))
            #sel.add(self.residue_range('212:E', '214:E'))

        if self.location == "preLC3":

            #sel = selection(self.residue_range('197:A', '199:A'))
            sel = selection(self.residue_range('196:C', '199:C'))
            #sel.add(self.residue_range('201:B', '203:B'))
            #sel.add(self.residue_range('201:D', '203:D'))
            #sel.add(self.residue_range('212:E', '214:E'))

        if self.location == "preLC4":

            #sel = selection(self.residue_range('197:A', '199:A'))
            #sel = selection(self.residue_range('197:C', '197:C'))
            #sel.add(self.residue_range('205:C', '207:C'))
            sel = selection(self.residue_range('202:B', '204:B'))
            sel.add(self.residue_range('202:D', '204:D'))
            #sel.add(self.residue_range('212:E', '214:E'))

        if self.location == "preLC1a":

            #sel = selection(self.residue_range('197:A', '199:A'))
            #sel.add(self.residue_range('197:C', '199:C'))
            #sel.add(self.residue_range('201:B', '203:B'))
            #sel.add(self.residue_range('201:D', '203:D'))
            sel = selection(self.residue_range('212:E', '213:E'))
            #sel.add(self.residue_range('220:E', '221:E'))

        if self.location == "preLC2a":

            sel = selection(self.residue_range('197:A', '199:A'))
            sel.add(self.residue_range('197:C', '199:C'))
            #sel.add(self.residue_range('201:B', '203:B'))
            #sel.add(self.residue_range('201:D', '203:D'))

        if self.location == "preLC2a_s":

            #sel = selection(self.residue_range('197:A', '199:A'))
            sel = selection(self.residue_range('197:C', '199:C'))
            sel.add(self.residue_range('205:C', '206:C'))
            #sel.add(self.residue_range('201:B', '203:B'))
            #sel.add(self.residue_range('201:D', '203:D'))

        if self.location == "preLC3a":

            sel = selection(self.residue_range('196:C', '197:C'))
            #sel.add(self.residue_range('205:C', '205:C'))

            #sel.add(self.residue_range('201:B', '203:B'))
            #sel.add(self.residue_range('201:D', '203:D'))

        return sel



