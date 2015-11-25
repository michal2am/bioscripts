set start_file "close_49_complete.pdb"
mol new $start_file

set all [atomselect top water]

puts "Calculating box dimensions"

set dim [measure minmax $all]
set dim0 [lindex $dim 0]
set dim1 [lindex $dim 1]

set xbox [expr { [lindex $dim1 0] -[lindex $dim0 0] } ]
set ybox [expr { [lindex $dim1 1] -[lindex $dim0 1] } ]
set zbox [expr { [lindex $dim1 2] -[lindex $dim0 2] } ]

set cent [measure center $all]
set xc [lindex $cent 0]
set yc [lindex $cent 1]
set zc [lindex $cent 2]

puts "Box dimensions: $xbox $ybox $zbox"
puts "Box center: $xc $yc $zc"

puts "Setting the box"

pbc set "$xbox $ybox $zbox"
pbc box -center origin
                                                     

