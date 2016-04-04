for {set i 0} {$i < 33} {incr i 3} {
[atomselect top protein frame $i] writepdb $i.pdb
}