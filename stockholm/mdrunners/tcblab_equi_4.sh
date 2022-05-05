#!/bin/bash

#SBATCH -p lindahl
#SBATCH -J 6x3z_wt_abu_equL
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -C gpu --gres=gpu:4
#SBATCH -e job-%j.err -o job-%j.out
#SBATCH -d singleton

rest_prefix=step5_input
mini_prefix=step6.0_minimization
equi_prefix=step6.%d_equilibration

cnt=1
cntmax=6

while (($cnt <= $cntmax))
do

        pcnt=$((cnt-1))

        istep=`printf ${equi_prefix} ${cnt}`
        pstep=`printf ${equi_prefix} ${pcnt}`

        if (($cnt==1)) ; then
                pstep=$mini_prefix
        fi

        module unload gromacs
        module load gromacs/2021.5

        sys=1
        sysmax=4

        while (($sys <= $sysmax))
        do
                cd sys${sys}
                gmx grompp -f ../mdp/${istep}.mdp -o ${istep}.tpr -c ${pstep}.gro -r ../${rest_prefix}.gro -p ../topol.top -n ../index.ndx
                cd ../
                sys=$((sys+1))
        done

        module unload gromacs
        module load gromacs/2021.5 gromacs=gmx_mpi

        mpirun -np 4 gmx_mpi mdrun -v -deffnm ${istep} -multidir sys1/ sys2/ sys3/ sys4/ -nb gpu -pme gpu -update gpu -bonded gpu

        cnt=$((cnt+1))

done
