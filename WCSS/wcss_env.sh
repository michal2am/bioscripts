#!/bin/bash
# Common environment for all WCSS jobs - sourced by wcss_mini.sh / wcss_equi.sh / wcss_prod.sh

module purge
module load GCC/13.3.0
module load cuda/12.8.0
module load OpenMPI/5.0.3-GCC-13.3.0
source $HOME/software/gromacs-2026.4/bin/GMXRC

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

GMX=gmx_mpi                       # MPI build is required for -multidir
#LAUNCH="srun --cpu-bind=none"     # GROMACS does the pinning (-pin on); alternative: "mpirun -np $SLURM_NTASKS"
LAUNCH="srun --mpi=pmix --cpu-bind=none" #trying to fix rank setting
MAXWARN=2                         # check which warnings these are, ideally go to 0

SYSTEMS=(sys1 sys2 sys3 sys4)     # each: step5_input.gro, topol.top, index.ndx, toppar/

command -v $GMX >/dev/null || { echo "ERROR: $GMX not found - check modules/GMXRC"; exit 1; }
module list 2>&1
