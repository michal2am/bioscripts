#!/bin/bash

#SBATCH -p lindahl
#SBATCH -J jobname
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -C gpu --gres=gpu:4
#SBATCH -e job-%j.err -o job-%j.out
#SBATCH -d singleton

init=step5_input
rest_prefix=step5_input
mini_prefix=step6.0_minimization

module unload gromacs
module load gromacs/2021.5

sys=1
sysmax=4

while (($sys <= $sysmax))
do
        cd sys${sys}/
        gmx grompp -f ../mdp/${mini_prefix}.mdp -o ${mini_prefix}.tpr -c ../${init}.gro -r ../${rest_prefix}.gro -p ../topol.top -n ../index.ndx
        cd ../
        sys=$((sys+1))
done

module unload gromacs
module load gromacs/2021.5 gromacs=gmx_mpi

mpirun -np 4 gmx_mpi mdrun -v -deffnm ${mini_prefix} -multidir sys1/ sys2/ sys3/ sys4/
