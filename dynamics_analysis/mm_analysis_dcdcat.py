from MDAnalysis import *

start = 5
stop = 15

psf = '06.gabar_ionized_2gaba.psf'
dcd = list()

fs2akma = 0.0204582651
ts = 1
efts = ts*fs2akma
savefreq = 25000

for i in range(start, stop+1):
        dcd.append('gabar_2gaba.%03d.dcd' % i)

print "Merging: "
print dcd

u = Universe(psf, dcd)
w = Writer('gabar_2gaba_all.dcd', u.atoms.numberOfAtoms(), delta=efts, step=savefreq )

for ts in u.trajectory:
        w.write(u)