

# Using readlines()
pdb = open('step1_pdbreader.pdb', 'r')
pdb_lines = pdb.readlines()


# Strips the newline character
for line in pdb_lines:
    if 'PROA' in line:
        print(line.strip().split())

    #print("Line{}: {}".format(count, line.strip()))


# writing to file
# file1 = open('myfile.txt', 'w')
# file1.writelines(L)
# file1.close()