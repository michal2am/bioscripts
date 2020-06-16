#!/bin/zsh

for f in *all_fit.pdb;
do
    echo "$f"
    python ~/repos/bioscripts/homology_modeling/model_evaluation/new_pl/model_distance.py -m $f -n 100 -p1 17 2896 -p2 8172 11051
done
