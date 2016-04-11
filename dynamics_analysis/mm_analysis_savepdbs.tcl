for {set i 0} {$i < 33} {incr i 3} {
[atomselect top protein frame $i] writepdb $i.pdb
}



for {set i 1000} {$i < 1010} {incr i 1} {
[atomselect top "name OH2 and (x*x + y*y < 50) and within 10 of pore" frame $i] get z
}



