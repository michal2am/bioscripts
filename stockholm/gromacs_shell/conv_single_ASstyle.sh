#!/bin/bash

# those should be done before to prepare files with selected atoms and chain IDs
# gmx make_ndx -f step5_input.gro -n index.ndx -o new_index.ndx
# gmx trajconv -f step5_input.gro -s md.tpr -n new_index.ndx -o new_pdb.pdb

pdbfile=$1    # new pdb file with selected atoms only
indexfile=$2  # new index file with selection corresponding to new pdb
select=$3     # new group number from new index
outname=$4    # should be the same as pdb and ndx, e.g. gaba_holo_protligpopcion_MD!.xtc

gmx trjconv -f traj_comp.xtc -s md.tpr -pbc whole -n $indexfile  -o view1.xtc <<EOF
$select
EOF

gmx trjconv -f view1.xtc -s $pdbfile -pbc nojump -n $indexfile -o view2.xtc <<EOF
$select
EOF

gmx trjconv -f view2.xtc -s md.tpr -pbc mol -ur compact -center -n $indexfile -o $outname <<EOF
3
$select
EOF

rm view1.xtc view2.xtc