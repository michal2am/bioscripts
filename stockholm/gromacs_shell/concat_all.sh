#! /bin/bash
# concatenates all converted .xtc files

gmx trjcat -f `ls step7_*.xtc | sort -n -t _ -k 2` -o all_runsteps.xtc -cat
