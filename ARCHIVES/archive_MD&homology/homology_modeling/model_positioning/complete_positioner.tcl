mol load pdb ../hybrid_models/hybrid_4cof_3rif/modeller_models_many/gabar_MM.01.B99990035_fit.pdb 
source orienter.tcl
source centerer.tcl
[atomselect top all] writepdb model_35_fit_pos.pdb
mol delete top
resetpsf