#! /bin/bash
# concatenates all raw xtc files from equilibration steps
# outname is for temporary file which has to be converted, eg. all_equisteps.xtc

outname=$1

gmx trjcat -f `ls step6.*.xtc | sort -n -t _ -k 2` -o "$outname" -cat
