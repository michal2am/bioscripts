#!/bin/bash

#SBATCH -p batch_16h
#SBATCH -J jobname
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 24
#SBATCH --output="job.out"
#SBATCH --error="job.err"

module load tryton/gromacs/2020.2

init=step5_input
rest_prefix=step5_input
mini_prefix=step6.0_minimization

sys=1
sysmax=4

while (($sys <= $sysmax))
do
        cd sys${sys}/
        gmx_mpi grompp -f ../mdp/${mini_prefix}.mdp -o ${mini_prefix}.tpr -c ../${init}.gro -r ../${rest_prefix}.gro -p ../topol.top -n ../index.ndx
        cd ../
        sys=$((sys+1))
done

mpiexec gmx_mpi mdrun -v -deffnm step6.0_minimization -multidir sys1/ sys2/ sys3/ sys4/