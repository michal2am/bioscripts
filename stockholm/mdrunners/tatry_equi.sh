#!/bin/bash

# set pinoffset  0, 8, 16, 24
# set gpu_id 0, 0, 1, 1
# nohup ./tatry_equi.sh > run.out 2> run.err &

init=step5_input
rest_prefix=step5_input
mini_prefix=step6.0_minimization
equi_prefix=step6.%d_equilibration

gmx grompp -f ../mdp/${mini_prefix}.mdp -o ${mini_prefix}.tpr -c ../${init}.gro -r ../${rest_prefix}.gro -p ../topol.top -n ../index.ndx
gmx mdrun -v -deffnm ${mini_prefix} -pin on -pinstride 1 -nt 8 -pinoffset 0 -gpu_id 0

cnt=1
cntmax=6

while (($cnt <= $cntmax))
do

        pcnt=$((cnt-1))

        istep=`printf ${equi_prefix} ${cnt}`
        pstep=`printf ${equi_prefix} ${pcnt}`

        if (($cnt==1)) ; then
                pstep=$mini_prefix
        fi

        gmx grompp -f ../mdp/${istep}.mdp -o ${istep}.tpr -c ${pstep}.gro -r ../${rest_prefix}.gro -p ../topol.top -n ../index.ndx -maxwarn 2
        gmx mdrun -v -deffnm ${istep} -nb gpu -pme gpu -update gpu -bonded gpu -pin on -pinstride 1 -ntmpi 1 -ntomp 8 -gpu_id 0 -pinoffset 0

        cnt=$((cnt+1))

done
