#! /bin/bash
# converts raw single .trr file from charmm-gui MD according to AS advices
# $1 is a file number according to standard charmm-gui notation and output is step7_$1.xtc
# test

step=$1

gmx trjconv -f step7_${step}.trr -s step7_${step}.tpr -pbc whole -o view1.xtc <<EOF
System
EOF

gmx trjconv -f view1.xtc -s step5_input.gro -pbc nojump -o view2.xtc <<EOF
System
EOF

gmx trjconv -f view2.xtc -s step7_${step}.tpr -pbc mol -ur compact -center -o step7_${step}.xtc <<EOF
Protein
System
EOF

rm view1.xtc view2.xtc