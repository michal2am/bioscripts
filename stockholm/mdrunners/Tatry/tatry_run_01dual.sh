#!/bin/bash

# nohup ./tatry_run_01dual.sh > run.out 2> run.err &

gmx grompp -f ../mdp/step7_production.mdp -o step7_production -c step6.6_equilibration.gro -p ../topol.top -n ../index.ndx
gmx mdrun -v -deffnm step7_production -cpi -pmefft gpu -nb gpu -pme gpu -update gpu -bonded gpu -pin on -pinstride 1 -ntmpi 1 -ntomp 16 -gpu_id 0 -pinoffset 0 -maxh 96
