#!/bin/bash
# prepares index and pdb file from charmm-gui files, just write down the selected group for trajectory work later on
# name template should be consistent, e.g 6x3z_gaba_protliglipion

name_template=$1
inname_tpr=$2

gmx make_ndx -f step5_input.gro -n index.ndx -o "$name_template".ndx
gmx trjconv -f step5_input.gro -s $inname_tpr -n "$name_template".ndx -o "$name_template".pdb