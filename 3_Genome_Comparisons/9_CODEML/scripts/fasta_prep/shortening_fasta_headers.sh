#!/bin/bash
#Description: Script for shortening and standardizing fasta headers

##################################### DEFINE FASTA AND NEW FASTA NAME #################################################
FASTA="Araport11_seq_20210622_representative_gene_model.fasta"
ABBR="Aratha" 

sed ':a;N;/^>/M!s/\n//;ta;P;D' $FASTA > UNWRAP_$FASTA
sed -i 's/\s.*$//' UNWRAP_$FASTA
sed -i 's/>/>'$ABBR'_/g' UNWRAP_$FASTA
mv UNWRAP_$FASTA $ABBR.SHORT.cds.fasta
