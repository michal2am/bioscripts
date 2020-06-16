#!/bin/zsh

for f in *log;
do
    echo "$f"
    python ~/repos/bioscripts/homology_modeling/model_evaluation/new_pl/model_eval_global.py -l $f  -n 100 -t 0.25
done