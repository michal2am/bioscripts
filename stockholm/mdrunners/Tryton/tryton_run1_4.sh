#!/bin/bash

#SBATCH -p batch_16h
#SBATCH -J jobname
#SBATCH -t 1:00:00
#SBATCH -N 4
#SBATCH --ntasks-per-node 24
#SBATCH --output="job.out"
#SBATCH --error="job.err"

module load tryton/gromacs/2020.2

sys=1
sysmax=4

while (($sys <= $sysmax))
do
        cd sys${sys}
        gmx_mpi grompp -f ../mdp/step7_production.mdp -o step7_production -c step6.6_equilibration.gro -p ../topol.top -n ../index.ndx
        cd ..
        sys=$((sys+1))
done

mpiexec gmx_mpi mdrun -v -deffnm step7_production -multidir sys1 sys2 sys3 sys4 -cpi -maxh 0.1