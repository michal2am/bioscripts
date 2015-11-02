#!/bin/bash

pdbDir="$1"
pdbAll="${pdbDir}/*.pdb"

echo "info: begin sequencing files in: ${pdbDir}"

for pdb in $pdbAll; do

        pdbFile="${pdb}"
        pdbName="$( basename "$pdbFile" | sed 's/\(.*\).pdb/\1/')"
        pdbPath="$( dirname "$pdbFile" )"

        newFastaFile="${pdbPath}/${pdbName}.fasta"

        echo "info: processing file "$pdbName" "

        python pdb2fasta.py --pdbFile "$pdbFile" --pdbCode "$pdbName" --seqName "$pdbName"

done

echo "info: sequencing done, moving files to: "${pdbDir}"/sequenced"

mkdir -p "$pdbDir"/sequenced
mv *.fasta "$pdbDir"/sequenced
