#!/bin/bash

#SBATCH -p lindahl
#SBATCH -J jobname
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -C gpu --gres=gpu:4
#SBATCH -e job-%j.err -o job-%j.out
#SBATCH -d singleton

module unload gromacs
module load gromacs/2021.5

sys=1
sysmax=4

while (($sys <= $sysmax))
do
        cd sys${sys}
        gmx grompp -f ../mdp/step7_production.mdp -o step7_production -c step6.6_equilibration.gro -p ../topol.top -n ../index.ndx
        cd ..
        sys=$((sys+1))
done

module unload gromacs
module load gromacs/2021.5 gromacs=gmx_mpi

mpirun -np 4 gmx_mpi mdrun -v -deffnm step7_production -multidir sys1 sys2 sys3 sys4 -cpi -maxh 0.1 -nb gpu -pme gpu -update gpu -bonded gpu