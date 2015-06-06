
class pdbLine:

    def __init__(self, atomLine):
        self.line = atomLine
        self.prefix = 'ATOM  '
        self.serial = int(atomLine[6:11])
        self.name = atomLine[12:16]
        self.resName = atomLine[17:20]
        self.chainID = atomLine[21]
        self.resSeq = int(atomLine[22:26])
        self.coors = [float(atomLine[30:38]), float(atomLine[38:46]), float(atomLine[46:54])]
        self.seg = atomLine[72:76]
	
    def writePDBline(self, serial=' ', chain=' ', res=' ', seg=' ', coors=' '):

        if serial == ' ':
            serial = self.serial
        if chain == ' ':		
            chain = self.chainID
        if res == ' ':
            res = self.resSeq
        if seg == ' ':
            seg = self.seg
        if coors == ' ':
           coors = [self.coors[i] for i in range(0,3)]

        newLine = '{0:s}{1:5d}{2:s} {3:s}{4:s}{5:4d}{6:s}{7:8.3f}{8:8.3f}{9:8.3f}{10:s}{11:s}\n'.format(self.prefix, serial, self.name, self.resName, chain, res, self.line[26:30], coors[0], coors[1], coors[2], self.line[54:72], seg)

        return newLine


class pdbFile:

    def __init__(self, old_pdbfile, new_pdbfile):
        old = open(old_pdbfile, 'r')
        new = open(new_pdbfile, 'w')
        old_lines = old.readlines()
        old_linestowrite = []
        for line in old_lines:
            if line[0:4] == 'ATOM':
                old_linestowrite.append(pdbLine(line).writePDBline())
        for line in old_linestowrite:
            new.write(line)




test = pdbFile('4COF_prot.pdb', 'new.pdb')
        

