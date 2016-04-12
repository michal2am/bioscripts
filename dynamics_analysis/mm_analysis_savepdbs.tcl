for {set i 0} {$i < 33} {incr i 3} {
[atomselect top protein frame $i] writepdb $i.pdb
}


for {set i 1} {$i < 15} {incr i 1} {
    set file [open "$i.dat" w]
    for {set j [expr {$i*100}]} {$j < [expr {$i*100 + 100}]} {incr j 1} {
    puts $file [[atomselect top "name OH2 and (x*x + y*y < 50) and within 10 of pore" frame $j] get z]
    }
    close $file
}





