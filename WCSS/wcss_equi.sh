#!/bin/bash -l
#SBATCH -J equi
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4                 # one MPI rank per system
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=4          # 1 = all 4 share a GPU; 4 = one GPU each (GROMACS assigns automatically)
#SBATCH --partition=lem-gpu
#SBATCH -o %x-%j.out

cd "$SLURM_SUBMIT_DIR" || exit 1
source ./wcss_env.sh
nvidia-smi

WALL_H=24          # keep in sync with -t above
MARGIN_H=0.25      # time left for writing checkpoints/outputs

mini=step6.0_minimization
equi=step6.%d_equilibration

all_done() { local s; for s in "${SYSTEMS[@]}"; do [[ -f $s/$1.gro ]] || return 1; done; }

for cnt in {1..6}; do
    istep=$(printf $equi $cnt)
    if (( cnt == 1 )); then pstep=$mini; else pstep=$(printf $equi $((cnt-1))); fi

    if all_done $istep; then echo "--- $istep already done ---"; continue; fi
    all_done $pstep || { echo "ERROR: $pstep.gro missing in some system"; exit 1; }

    # Exact continuation from previous step's checkpoint (minimization has no .cpt)
    tflag=""
    (( cnt > 1 )) && tflag="-t $pstep.cpt"

    for s in "${SYSTEMS[@]}"; do
        [[ -f $s/$istep.tpr ]] && continue
        ( cd "$s" && $GMX grompp -f ../mdp/$istep.mdp \
                                 -o $istep.tpr \
                                 -c $pstep.gro $tflag \
                                 -r step5_input.gro \
                                 -p topol.top \
                                 -n index.ndx \
                                 -maxwarn $MAXWARN ) || { echo "ERROR: grompp $istep failed in $s"; exit 1; }
    done

    # -maxh counts from mdrun start, so give each step only the time left in this job
    maxh=$(awk -v w=$WALL_H -v m=$MARGIN_H -v e=$SECONDS 'BEGIN{printf "%.2f", w-m-e/3600}')
    awk -v h=$maxh 'BEGIN{exit !(h>0.1)}' || { echo "Out of wall time before $istep - resubmit"; exit 1; }

    echo "--- Running $istep (maxh $maxh) ---"
    $LAUNCH $GMX mdrun -multidir "${SYSTEMS[@]}" \
        -deffnm $istep \
        -cpi \
        -nb gpu -pme gpu -bonded gpu -update gpu \
        -ntomp $OMP_NUM_THREADS \
        -pin on \
        -maxh $maxh

    all_done $istep || { echo "$istep not finished (maxh or crash) - resubmit this script, it continues from .cpt"; exit 1; }
done

echo "Equilibration done for all systems."
