#!/bin/bash
Fasta_File=$1
Project_name=$2
GFF3_File=$3

Name_size=$(echo $(head -n 1 $Fasta_File.old | wc -m )-2 |bc)

echo $Fasta_File
echo $Name_size
echo $Project_name
echo $GFF3_File

if [ ! -f "$GFF3_File" ]; then
   echo "$GFF3_File does not exist.....exiting with error."
   exit
fi

if [[ $Name_size -ge 15 ]]
then
  echo "The Name of the fasta file is greater than 15: Changing the name of GFF3 to match FASTA"
  mv $Fasta_File $Fasta_File.replace
  mv $Fasta_File.old $Fasta_File

  awk -v Project_name_awk="${Project_name}_Chr" '{gsub("_Chr",Project_name_awk);print}' $GFF3_File > $GFF3_File.corrected
  mv $GFF3_File $GFF3_File.replace
  mv $GFF3_File.corrected $GFF3_File
fi
