#!/bin/bash
# execute in a system directory (e.g. 6x3x_wt_ABU_DZP_equilong) containing charmm-gui system .tar.gz archive
# suitable mdp files should be copied manually to the mdp directory

system_archive=$1

mkdir full_system
tar xfzv $system_archive --strip-components=1 -C full_system/

cp -r full_system/gromacs/toppar/ .
cp full_system/gromacs/*.top .
cp full_system/gromacs/*.ndx .
cp full_system/gromacs/step5* .

# CG only:
cp full_system/gromacs/*.itp .

mkdir mdp
mkdir sys1 sys2 sys3 sys4