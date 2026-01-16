#!/bin/bash
#miniscript for making a directory with index files from Orthogroups.txt
#requires Orthogroups.txt in PARENT
#Define paths to PARENT directory
PARENT=/usit/abel/u1/siribi/nobackup/COLD/PIPELINE_RUN_DECEMBER2018/Building_alignments_for_CODEML/draba_MCOGs/

#Make index files in new directory called index_files
cd $PARENT
sed 's/://' Orthogroups.txt > Orthogroups_noColon.txt
grep "DnA" Orthogroups_noColon.txt | grep "DnN" | grep "Athal" | grep "Avern" | grep "Aalp" | grep "Acan" | grep "Dhis" | grep "Dnem" > Orthogroups_atleastoneofeach.txt
mkdir index_files
cd index_files
cp ../Orthogroups_atleastoneofeach.txt .
awk '{for(i=2;i<=NF;i++){printf "%s ", $i >> $1".txt"};printf "\n" >> $1".txt"; close($1".txt")}' Orthogroups_atleastoneofeach.txt
rm Orthogroups_atleastoneofeach.txt 
echo done! 
