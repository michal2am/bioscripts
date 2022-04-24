#!/bin/bash

# set pinoffset  0, 8, 16, 24
# set gpu_id 0, 0, 1, 1
# nohup ./tatry_run.sh > run.out 2> run.err &

gmx grompp -f ../mdp/step7_production.mdp -o step7_production -c step6.6_equilibration.gro -p ../topol.top -n ../index.ndx
gmx mdrun -v -deffnm step7_production -cpi -nb gpu -pme gpu -update gpu -bonded gpu -pin on -pinstride 1 -ntmpi 1 -ntomp 8 -gpu_id 0 -pinoffset 0 -maxh 6
