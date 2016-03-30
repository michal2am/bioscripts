for {set i 0} {$i < 2201} {incr i 200} {
[atomselect top protein frame $i] writepdb $i.pdb
}