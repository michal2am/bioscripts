# parameters
set inPdb "close_49_autopsf_pos.pdb"
set inPsf "close_49_autopsf.psf"
set outPdb "close_49_autopsf_pos_sys.pdb"
set outPsf "close_49_autopsf_pos_sys.psf"


# set echo on for debugging
echo on

# need psfgen module and topology
package require psfgen
topology top_all27_prot_lipid_na.inp

# load structures
resetpsf
readpsf $inPsf
coordpdb $inPdb
readpsf poun_x4.psf
coordpdb poun_x4.minim.coor.pdb


# write temporary structure
set temp "temp"
writepsf $temp.psf
writepdb $temp.pdb

# reload full structure (do NOT resetpsf!)
mol load psf $temp.psf pdb $temp.pdb

# select and delete waters that overlap protein:
set selwat [atomselect top "resname TIP3"]
set lseglist [lsort -unique [$selwat get segid]]
foreach lseg $lseglist {
 set selover [atomselect top "segid $lseg and within 5.1 of protein"]
 set resover [lsort -unique [$selover get resid]]
 foreach res $resover {
  delatom $lseg $res
 }
}
foreach res { } {delatom $WAT1 $res}
foreach res { } {delatom $WAT2 $res}

# select and delete lipids that overlap protein:
set selwat [atomselect top "resname POUN"]
set lseglist [lsort -unique [$selwat get segid]]
foreach lseg $lseglist {
 set selover [atomselect top "segid $lseg and within 2.2 of protein"]
 set resover [lsort -unique [$selover get resid]]
 foreach res $resover {
  delatom $lseg $res
 }
}
foreach res { } {delatom $WAT1 $res}
foreach res { } {delatom $WAT2 $res}

# write full structure
writepsf $outPsf
writepdb $outPdb

# clean up
file delete $temp.psf
file delete $temp.pdb

# non-interactive script:?? vmd -dispdev text < combine.tcl > combine.log
quit
 
