#!/bin/bash

#SBATCH -J equi
#SBATCH -t 24:00:00
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


rest_prefix=step5_input
mini_prefix=step6.0_minimization
equi_prefix=step6.%d_equilibration

for ((cnt=1; cnt<=6; cnt++)); do

    pcnt=$((cnt-1))
    istep=$(printf ${equi_prefix} ${cnt})
    pstep=$(printf ${equi_prefix} ${pcnt})

    if [ $cnt -eq 1 ]; then
        pstep=$mini_prefix
    fi

    if [ -f "sys4/${istep}.gro" ]; then
        echo "--- Step ${istep} already done. Starting next... ---"
        continue
    fi

    echo "--- Starting: $istep (previous: $pstep) ---"

    for sys in {1..3}; do
        cd sys${sys}
        if [ ! -f "${istep}.tpr" ]; then
            echo "--- Starting $istep with new .tpr file"
            gmx_mpi grompp -f ../mdp/${istep}.mdp \
                           -o ${istep}.tpr \
                           -c ${pstep}.gro \
                           -r ../${rest_prefix}.gro \
                           -p ../topol.top \
                           -n ../index.ndx \
                           -maxwarn 2
        fi
        cd ../
    done

    mpirun gmx_mpi mdrun -v \
        -multidir sys1/ sys2/ sys3/ \
        -deffnm ${istep} \
        -cpi \
        -nb gpu \
        -pme gpu \
        -pmefft gpu \
        -update gpu \
        -bonded gpu \
        -ntomp $OMP_NUM_THREADS \
        -pin on \
        -gpu_id 0 \
        -maxh 23.5
done
