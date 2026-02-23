
from MDAnalysis.analysis import hole2
import MDAnalysis as mda
import matplotlib.pyplot as plt



test = mda.Universe('step5_charmm2namd.psf', 'f65g_run_all_s10.dcd')
ha = hole2.HoleAnalysis(test, executable='~/hole2/exe/hole')
ha.run()

print(ha)
print(ha.profiles)
ha.plot()
plt.show()