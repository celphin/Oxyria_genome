#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=0-1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=120000m

#to run sbatch gff3_to_cds.sh gff3name.gff3 fastaname.fasta outputname.fasta

module load bedtools/2.30.0
module load bedops/2.4.41

gff3_file=$1
full_fasta=$2
output_name=$3

#Converts to bed
srun gff2bed < $gff3_file > annotation.bed
#Extracts CDS from the gff3 files
grep "gene" annotation.bed > cleaned_annotation.bed
#Get fasta based on bed file
srun bedtools getfasta -fi $full_fasta -bed cleaned_annotation.bed -fo mapped_fasta.fasta
awk '
  BEGIN { RS=">"; ORS="" }
  NR>1 { seq=$0; gsub("\n", "", seq); if (!sequences[seq]++) print ">"$0 }
' mapped_fasta.fasta >  $output_name
rm annotation.bed
rm cleaned_annotation.bed
rm mapped_fasta.fasta
