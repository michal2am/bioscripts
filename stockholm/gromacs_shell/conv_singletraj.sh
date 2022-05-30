#!/bin/bash
# ~/repos/bioscripts/stockholm/gromacs_shell/conv_singletraj.sh ../6x3z_eqXL.pdb ../6x3z_eqXL.ndx 3 5 step7_production.xtc step7_production.tpr 6x3z_eqXL_MD1.xtc

pdbfile=$1    # new pdb file with selected atoms only
indexfile=$2  # new index file with selection corresponding to new pdb
select=$3     # output group
center=$4     # center group
inname_xtc=$5 # step7_production.xtc
inname_tpr=$6 # step7_production.tpr
outname=$7    # should be the same as pdb and ndx, e.g. gaba_holo_protligpopcion_MD1.xtc

gmx trjconv -f $inname_xtc -s $inname_tpr -pbc whole -n $indexfile -o view1.xtc <<EOF
$select
EOF

# this one returns perfectly centered, but lipids diffused outside the box
gmx trjconv -f view1.xtc -s $pdbfile -pbc nojump -n $indexfile -o view2.xtc <<EOF
$select
EOF

# change -s tpr for pdb?
gmx trjconv -f view2.xtc -s $inname_tpr -pbc mol -ur compact -center -n $indexfile -o $outname <<EOF
$center
$select
EOF

rm view1.xtc view2.xtc


#gmx trjconv -f view2.xtc -s $inname_tpr -pbc mol -ur compact -center -n $indexfile -o $outname <<EOF
#3
#$select
#EOF