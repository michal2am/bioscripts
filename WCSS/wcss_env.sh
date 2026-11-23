#!/bin/bash
# Common environment for all WCSS jobs - sourced by wcss_mini.sh / wcss_equi.sh / wcss_prod.sh
# TODO: fill in module names once you have WCSS access (check: module avail gcc cuda mpi gromacs)

module purge
# module load gcc/XXX
# module load cuda/XXX
# module load <intel-oneapi-mpi | openmpi>/XXX
# Either your own build:
# source $HOME/software/gromacs_2026_gpu/bin/GMXRC
# ...or a site module:
# module load gromacs/XXX

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

GMX=gmx_mpi                       # MPI build is required for -multidir
LAUNCH="srun --cpu-bind=none"     # GROMACS does the pinning (-pin on); alternative: "mpirun -np $SLURM_NTASKS"
MAXWARN=2                         # check which warnings these are, ideally go to 0

SYSTEMS=(sys1 sys2 sys3 sys4)     # each: step5_input.gro, topol.top, index.ndx, toppar/

command -v $GMX >/dev/null || { echo "ERROR: $GMX not found - check modules/GMXRC"; exit 1; }
module list 2>&1
