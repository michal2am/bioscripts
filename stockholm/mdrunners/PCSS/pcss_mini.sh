#!/bin/bash

#SBATCH -J mini
#SBATCH -t 2:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 3
#SBATCH --cpus-per-task 10

module purge
module load gcc/14.2.0-gcc-11.5.0
module load cuda/12.8.0_570.86.10
module load intel-oneapi-mpi/2021.14.0-gcc-14.2.0
source $HOME/software/gromacs_2026_gpu/bin/GMXRC
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module list

init=step5_input
rest_prefix=step5_input
mini_prefix=step6.0_minimization

for sys in {1..3}; do
    cd sys${sys}/
    gmx_mpi grompp -f ../mdp/${mini_prefix}.mdp \
                   -o ${mini_prefix}.tpr \
                   -c ../${init}.gro \
                   -r ../${rest_prefix}.gro \
                   -p ../topol.top \
                   -n ../index.ndx \
                   -maxwarn 2
    cd ../
done

mpirun -np 3 gmx_mpi mdrun -v -deffnm ${mini_prefix} -multidir sys1/ sys2/ sys3/ -ntomp $OMP_NUM_THREADS -pin on