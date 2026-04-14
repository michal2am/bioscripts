## "E153_R207_BS1", "E153_R207_BS2", "E155_R207_BS1", "E155_R207_BS2"

# R207 site
# E153-R207 (R207 site)
set sel1 [atomselect top "chain A and resid 153 and name CD"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain A and resid 207 and name NH1"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

set sel1 [atomselect top "chain C and resid 153 and name CD"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain C and resid 207 and name NH1"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

# E155-R207 (R207 site)
set sel1 [atomselect top "chain A and resid 155 and name CD"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain A and resid 207 and name NH2"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

set sel1 [atomselect top "chain C and resid 155 and name CD"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain C and resid 207 and name NH2"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete