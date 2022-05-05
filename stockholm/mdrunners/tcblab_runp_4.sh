#!/bin/bash

#SBATCH -p lindahl
#SBATCH -J 6x3z_wt_abu
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -C gpu --gres=gpu:4
#SBATCH -e job-%j.err -o job-%j.out
#SBATCH -d singleton

module unload gromacs
module load gromacs/2021.5 gromacs=gmx_mpi

mpirun -np 4 gmx_mpi mdrun -v -deffnm step7_production -multidir sys1 sys2 sys3 sys4 -cpi -maxh 23.5 -nb gpu -pme gpu -update gpu -bonded gpu
