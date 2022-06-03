#!/bin/bash
# prepares index and pdb file from charmm-gui files, just write down the selected group for trajectory work later on
# name template should be consistent, e.g 6x3z_gaba_protliglipion

name_template=$1
inname_str=$2
inname_tpr=$3

gmx make_ndx -f $inname_str -n index.ndx -o "$name_template".ndx
gmx trjconv -f $inname_str -s $inname_tpr -n "$name_template".ndx -o "$name_template".pdb