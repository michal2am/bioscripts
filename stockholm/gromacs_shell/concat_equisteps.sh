#! /bin/bash
# concatenates all raw xtc files from equilibration steps, timestamps are cumulative durations + cumulative 50 ps (xtc step size)
# concat_equisteps.sh 8si9_ABUY4B_equi_sys1_newEqui.xtc

outname=$1


cp step6.1_equilibration.xtc step6.1_equilibration_timestamp.xtc
gmx trjconv -f step6.2_equilibration.xtc -o step6.2_equilibration_timestamp.xtc -t0 2550   # in ps
gmx trjconv -f step6.3_equilibration.xtc -o step6.3_equilibration_timestamp.xtc -t0 5100   # in ps
gmx trjconv -f step6.4_equilibration.xtc -o step6.4_equilibration_timestamp.xtc -t0 10150   # in ps
gmx trjconv -f step6.5_equilibration.xtc -o step6.5_equilibration_timestamp.xtc -t0 20200   # in ps
gmx trjconv -f step6.6_equilibration.xtc -o step6.6_equilibration_timestamp.xtc -t0 40250   # in ps

#1000350

gmx trjcat -f `ls step6.*_timestamp.xtc | sort -n -t _ -k 2` -o "$outname" -cat
