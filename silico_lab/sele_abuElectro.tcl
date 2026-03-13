# ABU carboxy
# T202
set sel1 [atomselect top "chain A and resid 202 and name OG1"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain H and resname ABU and name C4"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

# R67
set sel1 [atomselect top "chain B and resid 67 and name CZ"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain H and resname ABU and name C4"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

# T130
set sel1 [atomselect top "chain B and resid 130 and name OG1"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain H and resname ABU and name C4"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete


# ABU amino
# E155
set sel1 [atomselect top "chain A and resid 155 and name CD"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain H and resname ABU and name N"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

# Y97
set sel1 [atomselect top "chain A and resid 97 and name OH"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain H and resname ABU and name N"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

# S156
set sel1 [atomselect top "chain A and resid 156 and name O"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain H and resname ABU and name N"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete