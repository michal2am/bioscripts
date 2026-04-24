#! /bin/bash
# concatenates already concatenated equisteps with production run, timestamps are cumulative durations + cumulative 50 ps (xtc step size)
# concat_equirunsteps.sh  8si9_ABUY4B_equi_sys1_newEqui.xtc 8si9_ABUY4B_all_sys1_newEqui.xtc

equisteps=$1
outname=$2

gmx trjconv -f step7_production.xtc -o step7_production_timestamp.xtc -t0 100300   # in ps
#1000350

gmx trjcat -f "$equisteps" step7_production_timestamp.xtc -o "$outname" -cat