#!/bin/bash
#miniscript for making a directory with index files from Orthogroups.txt
#requires Orthogroups.txt in PARENT
#Define paths to PARENT directory
PARENT=/cluster/home/siribi/work/EVOTREE/CODEML/alignments/All40/

#Make index files in new directory called index_files
cd $PARENT
grep "Poptre" Orthogroups.tsv | grep "Popides" | grep "Salpur" | grep "Salsuc" | grep "Avecar" | grep "Cansat" | grep "Ficmic" | grep "Distri" | grep "Erijap" | grep "Malbac" | grep "Prunavi" | grep "Prunper" | grep "Pyrbet" | grep "Cucsat" | grep "Cucarg" | grep "Momcha" | grep "Ampedge" | grep "Cercan" | grep "Chafas" | grep "Dalodo" | grep "Faialb" | grep "Glymax" | grep "Phavul" | grep "Vigung" | grep "Alnglu" | grep "Betpen" | grep "Carfan" | grep "Casgla" | grep "Casmol" | grep "Fagsyl" | grep "Querob" | grep "Carill" | grep "Jugreg" | grep "Jatcur" | grep "Merann" | grep "Linusi" | grep "Kanobo" | grep "Aratha" | grep "Vitvin" | grep "Macint" > Orthogroups_atleastoneofeach.txt
mkdir index_files
cd index_files
cp ../Orthogroups_atleastoneofeach.txt .
awk '{for(i=2;i<=NF;i++){printf "%s ", $i >> $1".txt"};printf "\n" >> $1".txt"; close($1".txt")}' Orthogroups_atleastoneofeach.txt
rm Orthogroups_atleastoneofeach.txt 
echo done! 
