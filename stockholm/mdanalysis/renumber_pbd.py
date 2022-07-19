# simple script to shift numbering in given chain by some offset starting from given residue
# made for messed up cif files with alternate numberings

pdb = open('step1_pdbreader.pdb', 'r')
pdb_lines = pdb.readlines()

chain = 'PROA'
last_correct = 209
shift = 368

new_pdb_lines = []
# Strips the newline character
for line in pdb_lines:
    if chain in line:
        if int(line.strip()[23:26]) > last_correct:
            print(line.strip())
            new_line = line.strip()[0:23] + str(int(line.strip()[23:26]) - shift) + line.strip()[26:] + '\n'
            print(new_line)
            new_pdb_lines.append(new_line)
        else:
            new_pdb_lines.append(line)
    else:
        new_pdb_lines.append(line)

new_pdb = open('step1_pdbreader_renum.pdb', 'w')
new_pdb.writelines(new_pdb_lines)
new_pdb.close()
