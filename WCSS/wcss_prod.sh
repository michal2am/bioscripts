#!/bin/bash -l
#SBATCH -J prod
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4                 # one MPI rank per system
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:hopper:4          # 1 = all 4 share a GPU; 4 = one GPU each
#SBATCH --mem=64G
#SBATCH --partition=lem-gpu
#SBATCH -o %x-%j.out
#
# Same script for the first run and every restart: grompp only if .tpr is missing,
# mdrun -cpi continues from the checkpoint. Chain with: sbatch -d afterany:<jobid> wcss_prod.sh

cd "$SLURM_SUBMIT_DIR" || exit 1
source ./wcss_env.sh
nvidia-smi

WALL_H=24          # keep in sync with -t above
MARGIN_H=0.25      # time left for writing checkpoints/outputs

prod=step7_production
prev=step6.6_equilibration

# Nothing to do if all systems finished
done_all=1
for s in "${SYSTEMS[@]}"; do [[ -f $s/$prod.gro ]] || done_all=0; done
(( done_all )) && { echo "Production finished in all systems."; exit 0; }

for s in "${SYSTEMS[@]}"; do
    [[ -f $s/$prod.tpr ]] && continue
    ( cd "$s" && $GMX grompp -f ../mdp/$prod.mdp \
                             -o $prod.tpr \
                             -c $prev.gro \
                             -t $prev.cpt \
                             -p topol.top \
                             -n index.ndx \
                             -maxwarn $MAXWARN ) || { echo "ERROR: grompp failed in $s"; exit 1; }
done

maxh=$(awk -v w=$WALL_H -v m=$MARGIN_H -v e=$SECONDS 'BEGIN{printf "%.2f", w-m-e/3600}')

$LAUNCH $GMX mdrun -multidir "${SYSTEMS[@]}" \
    -deffnm $prod \
    -cpi \
    -nb gpu -pme gpu -bonded gpu -update gpu \
    -ntomp $OMP_NUM_THREADS \
    -pin on \
    -maxh $maxh
