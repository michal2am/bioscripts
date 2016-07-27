import mm_lib_filer as mmfil
import subprocess as sp

where = '.'
sext = '.FOR'
oext = '.o'
main = './SCAN.FOR'

sources = mmfil.find_files(where, sext)

for source in sources:

    if main == source:
        print("Omitting main file: {}".format(source))

    else:
        print("Compiling: {}".format(source))
        sp.run(['gfortran', '-c', '-std=f95', source])

objects = mmfil.find_files(where, oext)
comp_all = ['gfortran', main]
comp_all.extend(objects)
comp_all.extend(['-o', 'test'])

print("Putting all together:")
print(comp_all)
sp.run(comp_all)
