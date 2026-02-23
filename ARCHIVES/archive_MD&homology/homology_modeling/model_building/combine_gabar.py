
extr_dom = {'A': '', 'B': '', 'C': '', 'D': '', 'E': ''}
intr_dom = {'A': '', 'B': '', 'C': '', 'D': '', 'E': ''}
#attach chains according to structural alignment!
#because chain name's of both parts do not fit after alignement
chains = [['A', 'D'], ['B', 'C'], ['C', 'B'], ['D', 'A'], ['E', 'E']]

for line in open('4COF3rif_final.pdb'):
	if line[0:4] == 'ATOM':
		extr_dom[line[21]] += line

for line in open('3RIF_final_short.pdb'):
        if line[0:4] == 'ATOM':
                intr_dom[line[21]] += line

newPdb = open ('chimeric.pdb', 'w')

for chain in chains:
	newPdb.write(extr_dom[chain[0]])
	newPdb.write(intr_dom[chain[1]])

newPdb.close()

