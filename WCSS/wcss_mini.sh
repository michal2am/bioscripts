#!/bin/bash -l
#SBATCH -J mini_test
#SBATCH -t 08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4                 # one MPI rank per system
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=1
##SBATCH --partition=lem-gpu       # on GPU because CUDA libs are needed for GPU-compiled gmx
#SBATCH -o %x-%j.out

cd "$SLURM_SUBMIT_DIR" || exit 1
source ./wcss_env.sh

mini=step6.0_minimization

for s in "${SYSTEMS[@]}"; do
    [[ -f $s/$mini.tpr ]] && continue
    ( cd "$s" && $GMX grompp -f ../mdp/$mini.mdp \
                             -o $mini.tpr \
                             -c step5_input.gro \
                             -r step5_input.gro \
                             -p topol.top \
                             -n index.ndx \
                             -maxwarn $MAXWARN ) || { echo "ERROR: grompp failed in $s"; exit 1; }
done

$LAUNCH $GMX mdrun -multidir "${SYSTEMS[@]}" \
    -deffnm $mini \
    -ntomp $OMP_NUM_THREADS \
    -pin on

for s in "${SYSTEMS[@]}"; do
    [[ -f $s/$mini.gro ]] || { echo "ERROR: minimization not finished in $s"; exit 1; }
done
echo "Minimization done for all systems."
