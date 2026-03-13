# B9-B10 strands h-bonds
# 207-196
set sel1 [atomselect top "chain A and resid 207 and name N"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain A and resid 196 and name O"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

set sel1 [atomselect top "chain A and resid 207 and name O"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain A and resid 196 and name N"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

# 205 - 198
set sel1 [atomselect top "chain A and resid 205 and name N"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain A and resid 198 and name O"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

set sel1 [atomselect top "chain A and resid 205 and name O"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain A and resid 198 and name N"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

# 203 - 200
set sel1 [atomselect top "chain A and resid 203 and name N"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain A and resid 200 and name O"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

set sel1 [atomselect top "chain A and resid 203 and name O"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain A and resid 200 and name N"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete