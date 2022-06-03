#!/bin/bash
# ~/repos/bioscripts/stockholm/gromacs_shell/conv_singletraj.sh ../6x3z_eqXL.pdb ../6x3z_eqXL.ndx 3 5 step7_production.xtc step7_production.tpr 6x3z_eqXL_MD1.xtc

pdbfile=$1    # new pdb file with selected atoms only
indexfile=$2  # new index file with selection corresponding to new pdb
select=$3     # output group
center=$4     # center group
inname_xtc=$5 # step7_production.xtc
inname_tpr=$6 # step7_production.tpr
outname=$7    # just the template, e.g. gaba_holo_protligpopcion_MD

cd sys1
~/repos/bioscripts/stockholm/gromacs_shell/conv_singletraj.sh ../$pdbfile ../$indexfile $select $center $inname_xtc $inname_tpr ${outname}1.xtc &
cd ../sys2
~/repos/bioscripts/stockholm/gromacs_shell/conv_singletraj.sh ../$pdbfile ../$indexfile $select $center $inname_xtc $inname_tpr ${outname}2.xtc &
cd ../sys3
~/repos/bioscripts/stockholm/gromacs_shell/conv_singletraj.sh ../$pdbfile ../$indexfile $select $center $inname_xtc $inname_tpr ${outname}3.xtc &
cd ../sys4
~/repos/bioscripts/stockholm/gromacs_shell/conv_singletraj.sh ../$pdbfile ../$indexfile $select $center $inname_xtc $inname_tpr ${outname}4.xtc &
cd ..
