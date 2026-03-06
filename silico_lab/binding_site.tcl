# ligands

mol representation CPK
mol selection {resname ABU}
mol color Element
mol material AOChalky
mol addrep top

# protein chains

mol representation NewCartoon
mol selection {protein and chain A}
mol color ColorID 31
mol material AOChalky
mol addrep top

mol representation NewCartoon
mol selection {protein and chain B}
mol color ColorID 23
mol material AOChalky
mol addrep top

mol representation NewCartoon
mol selection {protein and chain C}
mol color ColorID 9
mol material AOChalky
mol addrep top

mol representation NewCartoon
mol selection {protein and chain D}
mol color ColorID 15
mol material AOChalky
mol addrep top

mol representation NewCartoon
mol selection {protein and chain E}
mol color ColorID 2
mol material AOChalky
mol addrep top

# protein residues

mol representation CPK
mol selection {chain B D and resid 46 67 and not element H}
mol color Element
mol material AOChalky
mol addrep top

mol representation CPK
mol selection {chain A C and resid 155 200 and not element H}
mol color Element
mol material AOChalky
mol addrep top

# bonds

set sel1 [atomselect top "chain A and resid 200 and name CG"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain B and resid 46 and name CG"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

set sel1 [atomselect top "chain C and resid 200 and name CG"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain D and resid 46 and name CG"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

set sel1 [atomselect top "chain A and resid 155 and name CD"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain H and resname ABU and name N"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

set sel1 [atomselect top "chain C and resid 155 and name CD"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain I and resname ABU and name N"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

set sel1 [atomselect top "chain B and resid 67 and name CZ"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain H and resname ABU and name C4"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

set sel1 [atomselect top "chain D and resid 67 and name CZ"]
set idx1 [$sel1 get index]
set sel2 [atomselect top "chain I and resname ABU and name C4"]
set idx2 [$sel2 get index]
label add Bonds 0/$idx1 0/$idx2
$sel1 delete
$sel2 delete

