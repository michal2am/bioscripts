#!/bin/bash

#SBATCH -p batch
#SBATCH -J jobname
#SBATCH -t 72:00:00
#SBATCH -N 4
#SBATCH --ntasks-per-node 24
#SBATCH --output="job.out"
#SBATCH --error="job.err"

module load tryton/gromacs/2020.2

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

        sys=1
        sysmax=4

        while (($sys <= $sysmax))
        do
                cd sys${sys}
                gmx_mpi grompp -f ../mdp/${istep}.mdp -o ${istep}.tpr -c ${pstep}.gro -r ../${rest_prefix}.gro -p ../topol.top -n ../index.ndx
                cd ../
                sys=$((sys+1))
        done

        mpiexec gmx_mpi mdrun -v -deffnm ${istep} -cpi -multidir sys1/ sys2/ sys3/ sys4/

        cnt=$((cnt+1))

done