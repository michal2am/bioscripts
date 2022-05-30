#!/bin/bash

#SBATCH -p lindahl
#SBATCH -J jobname
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -C gpu --gres=gpu:4
#SBATCH -e job-%j.err -o job-%j.out
#SBATCH -d singleton

sys=1
sysmax=4

# minimization step 0
# needs to be done using double precision
# TODO: make it parallel

while (($sys <= $sysmax))
do
        cd sys${sys}/

        module unload gromacs
        module load gromacs/2021.5

        gmx grompp -f ../mdp/step6.0_minimization.mdp -o step6.0_minimization.tpr -c ../step5_charmm2gmx.pdb -r ../step5_charmm2gmx.pdb -p ../system.top -n ../index.ndx -maxwarn -1

        module unload gromacs
        module load gromacs/2021.5 gromacs=gmx_d

        gmx_d mdrun -v -deffnm step6.0_minimization

        cd ../
        sys=$((sys+1))
done

# minimization step 1

sys=1

module unload gromacs
module load gromacs/2021.5

while (($sys <= $sysmax))
do
        cd sys${sys}/
        gmx grompp -f ../mdp/step6.1_minimization.mdp -o step6.1_minimization.tpr -c step6.0_minimization.gro -r ../step5_charmm2gmx.pdb -p ../system.top -n ../index.ndx -maxwarn -1
        cd ../
        sys=$((sys+1))
done

module unload gromacs
module load gromacs/2021.5 gromacs=gmx_mpi

mpirun -np 4 gmx_mpi mdrun -v -deffnm step6.1_minimization -multidir sys1/ sys2/ sys3/ sys4/