#!/bin/bash

#SBATCH -J run
#SBATCH -t 0:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=10
#SBATCH --gpus-per-node=1
#SBATCH --partition=tesla

module purge
module load gcc/14.2.0-gcc-11.5.0
module load cuda/12.8.0_570.86.10
module load intel-oneapi-mpi/2021.14.0-gcc-14.2.0
source $HOME/software/gromacs_2026_gpu/bin/GMXRC
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module list
nvidia-smi

for sys in {1..4}; do
    cd sys${sys}
    gmx_mpi grompp -f ../mdp/step7_production.mdp \
                   -o step7_production \
                   -c step6.6_equilibration.gro \
                   -p ../topol.top \
                   -n ../index.ndx \
                   -maxwarn 2
    cd ..
done

mpirun gmx_mpi mdrun -v \
    -multidir sys1/ sys2/ sys3/ \
    -deffnm step7_production \
    -cpi \
    -nb gpu \
    -pme gpu \
    -pmefft gpu \
    -update gpu \
    -bonded gpu \
    -ntomp $OMP_NUM_THREADS \
    -pin on \
    -gpu_id 0 \
    -maxh 0.1
