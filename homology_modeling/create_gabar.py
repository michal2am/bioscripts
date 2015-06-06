from modeller import *
from modeller.automodel import *   
from modeller.parallel import *
from gaba_modeller import *


j = job()
j.append(local_slave())
j.append(local_slave())
j.append(local_slave())
j.append(local_slave())

log.verbose()   
env = environ() 
env.io.atom_files_directory = ['.']

a = GabaModel(env,
              alnfile  = 'gabar_MM.01chimeric.pir',
              knowns   = 'chimeric',    
              sequence = 'gabar_MM.01',
	      assess_methods = (assess.DOPE, assess.GA341))            

a.starting_model= 1                 
a.ending_model  = 64
a.initial_malign3d = False
a.final_malign3d = True
a.md_level = refine.very_slow
a.library_schedule = autosched.slow  
a.use_parallel_job(j)                 
a.make()                            
a.cluster
