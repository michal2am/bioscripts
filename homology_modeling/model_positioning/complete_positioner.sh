#!/bin/bash

orienter='/home/mm/PycharmProjects17/bioscripts/homology_modeling/model_positioning/orienter.tcl'
centerer='/home/mm/PycharmProjects17/bioscripts/homology_modeling/model_positioning/centerer.tcl'

pdbDir="$1"
pdbAll="${pdbDir}/*.pdb"

echo "info: begin positioning files in: ${pdbDir}"

for pdb in $pdbAll; do

	pdbFile="${pdb}"
        pdbName="$( basename "$pdbFile" | sed 's/\(.*\).pdb/\1/')"
        pdbPath="$( dirname "$pdbFile" )"
 
	newPdbFile="${pdbPath}/${pdbName}_pos.pdb"

        echo "info: processing file "$pdbName" "


	vmdConfig="mol load pdb "$pdbFile"\n\
source "$orienter"\n\
source "$centerer"\n\
[atomselect top protein] writepdb "$newPdbFile"\n\
mol delete top\n\
resetpsf\n\
exit"

	echo -e "$vmdConfig" > vmdConfig.tcl
        
        vmd -e vmdConfig.tcl

done

echo "info: positioning done, moving files to: "${pdbDir}"/positioned"

mkdir -p "$pdbDir"/positioned
mv "*pos.pdb" "$pdbDir"/positioned
