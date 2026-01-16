#!/bin/bash
#miniscript for making a directory with index files from Orthogroups.txt
#requires Orthogroups.txt in PARENT
#Define paths to PARENT directory
PARENT=/cluster/home/siribi/work/EVOTREE/CODEML/alignments/MALPIGH/

#Make index files in new directory called index_files
cd $PARENT
grep "Merann" Orthogroups.tsv | grep "Jatcur" | grep "Kanobo" | grep "Salsuc" | grep "Salpur" | grep "Popides" | grep "Poptre" | grep "Linusi" | grep "Aratha" | grep "Vitvin" > Orthogroups_atleastoneofeach.txt
mkdir index_files
cd index_files
cp ../Orthogroups_atleastoneofeach.txt .
awk '{for(i=2;i<=NF;i++){printf "%s ", $i >> $1".txt"};printf "\n" >> $1".txt"; close($1".txt")}' Orthogroups_atleastoneofeach.txt
rm Orthogroups_atleastoneofeach.txt 
echo done! 
