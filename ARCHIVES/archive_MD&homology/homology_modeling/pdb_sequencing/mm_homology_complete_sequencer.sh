#!/bin/bash
# gets sequences of all pdb files in given directory
# michaladammichalowski@gmail.com
# ? - creation
# 14.12.15 - refactor
# EXAMPLE CALL: ./mm_homology_complete_sequencer.sh some_dir

pdbDir="$1"
pdbAll="${pdbDir}/*.pdb"
converter="${HOME}/pycharm_projects/bioscripts/homology_modeling/model_alignment/mm_homology_pdb2fasta.py"

echo "info: begin sequencing files in: ${pdbDir}"

for pdb in $pdbAll; do

        pdbFile="${pdb}"
        pdbName="$( basename "$pdbFile" | sed 's/\(.*\).pdb/\1/')"
        pdbPath="$( dirname "$pdbFile" )"

        newFastaFile="${pdbPath}/${pdbName}.fasta"

        echo "info: processing file "$pdbName" "

        python "$converter" --pdbFile "$pdbFile" --pdbCode "$pdbName" --seqName "$pdbName"

done

echo "info: sequencing done, moving files to: "${pdbDir}"/sequenced"

mkdir -p "$pdbDir"/sequenced
mv *.fasta "$pdbDir"/sequenced
