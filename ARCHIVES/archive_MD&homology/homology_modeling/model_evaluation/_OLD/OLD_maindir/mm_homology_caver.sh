#!/bin/bash

# run script for caver 3.0
# michaladammichalowski@gmail.com
# 26.12.15 - creation
# EXAMPLE CALL: ./mm_homology_caver.sh --pdbs . --conf mm_homology_conf.txt --outdir mod56 --caverdir /media/mm/data/soft/caver_3.0


while [[ $# > 1 ]]
do
key="$1"

case $key in
    -p|--pdbs)
    pdbs="$2"
    shift
    ;;
    -c|--conf)
    conf="$2"
    shift
    ;;
    -o|--outdir)
    outd="$2"
    shift
    ;;
    -c|--caverdir)
    cadr="$2"
    shift
    ;;
    *)
    echo "Not knonwn argument"
    ;;
esac
shift
done

caver_lib="${cadr}/caver/lib"
caver_jar="${cadr}/caver/caver.jar"
caver_hom="${cadr}/caver"

java -Xmx4000m         \
     -cp $caver_lib    \
     -jar $caver_jar   \
     -home $caver_hom  \
     -pdb $pdbs        \
     -conf $conf       \
     -out $outd