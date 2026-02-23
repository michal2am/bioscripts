set all [atomselect top protein]
set old [measure center $all]

set xmov [expr {-[lindex $old 0]} ]
set ymov [expr {-[lindex $old 1]} ]
set zmov [expr {-[lindex $old 2]} ]

set amov "$xmov $ymov $zmov"

$all moveby $amov

