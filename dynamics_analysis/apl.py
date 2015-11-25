#script for calculating area per lipid and system thickness
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')


xst = open('minim_patch32.001.xst', 'r').readlines()
xst.extend(open('minim_patch32.002.xst', 'r').readlines()[2:])
xst.extend(open('minim_patch32.003.xst', 'r').readlines()[2:])


leafSize = 32
time = 0.0
dtime = 2.0*1000/1000000	#timestep*framesave/fs2ns

timeAll = []
aplAll = []
thcAll = []
symAll = []

for line in xst[2:]:
	line = line.split()
	apl = float(line[1])*float(line[5])/leafSize
	sym = float(line[1])/float(line[5])
	thc = float(line[9])
#	print '%d %f %f' % (int(line[0]), apl, thc)
	
	timeAll.append(time)
	time = time + dtime
	aplAll.append(apl)
	thcAll.append(thc)
	symAll.append(sym)

fig = plt.figure(figsize=(6, 8))
ax = fig.add_subplot(311)
ax.plot(timeAll, aplAll, 'k-',  label="area per lipid")
ax.legend(loc="best")
ax.set_ylabel("APL [A]")
ax.set_xlabel("time [ns]")
ay = fig.add_subplot(312)
ay.plot(timeAll, thcAll, 'k--', label="system thickness")
ay.legend(loc="best")
ay.set_ylabel("thickness [A]")
ay.set_xlabel("time [ns]")
az = fig.add_subplot(313)
az.plot(timeAll, symAll, 'k--', label="system symmetry")
az.legend(loc="best")
az.set_ylabel("symmetry X:Y")
az.set_xlabel("time (ps)")

fig.savefig("apl.pdf")
