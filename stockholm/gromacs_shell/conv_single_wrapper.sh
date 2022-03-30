#! /bin/bash
# wrapper for 'conv_single.sh', $1 and $2 are run numbers according to standard charmm-gui notation

step=$1
laststep=$2

while [ $step -lt $laststep ]
do
        echo $step
        ./conv_single.sh $step
        step=$((step+1))
done