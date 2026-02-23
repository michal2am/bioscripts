# set pinoffset  0, 8, 16, 24
# set gpu_id 0, 0, 1, 1
# nohup ./local_equi.sh > run.out 2> run.err &

init=step5_input
rest_prefix=step5_input
mini_prefix=step6.0_minimization
equi_prefix=step6.%d_equilibration

gmx grompp -f ../mdp/${mini_prefix}.mdp -o ${mini_prefix}.tpr -c ../${init}.gro -r ../${rest_prefix}.gro -p ../topol.top -n ../index.ndx
gmx mdrun -v -deffnm ${mini_prefix} -pin on -pinstride 1 -nt 10 -pinoffset 0