#!/bin/bash

#SBATCH -p batch
#SBATCH -J jobname
#SBATCH -t 72:00:00
#SBATCH -N 4
#SBATCH --ntasks-per-node 24
#SBATCH --output="job.out"
#SBATCH --error="job.err"

module load tryton/gromacs/2020.2

mpiexec gmx_mpi mdrun -v -deffnm step7_production -multidir sys1 sys2 sys3 sys4 -cpi -maxh 72.5