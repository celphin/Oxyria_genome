#!/bin/bash
Fasta_File=$1
Project_name=$2

Name_size=$(echo $(head -n 1 $Fasta_File | wc -m )-2 |bc)

echo $Fasta_File
echo $Name_size
echo $Project_name
if [[ $Name_size -ge 15 ]]
then
  echo "The Name of the fasta file is greater than 15: Correcting name for EDTA"
  awk -v Project_name_awk="$Project_name" '{gsub(Project_name_awk,"");print}' $Fasta_File > $Fasta_File.replace 
  mv $Fasta_File $Fasta_File.old
  mv $Fasta_File.replace $Fasta_File
fi
