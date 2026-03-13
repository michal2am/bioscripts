
# loop C cap
# T202-Y205 h-bond
set sel1 [atomselect top "chain A and resid 205 and name OH"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain A and resid 202 and name OG1"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

# F200 - F46
set sel1 [atomselect top "chain A and resid 200 and name CG"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain B and resid 46 and name CG"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete