#!/bin/bash
#Script to prep files for mafft alignments
# Usage: ./run_pal2nal.sh 
# execute in directory with fasta files
module load EMBOSS/6.6.0-foss-2018b
FALIST=`ls -1 *.fasta`
i=1
for FAFILE in $FALIST
        do
	FILEBASE=$(echo $FAFILE | cut -f1 -d.)
	echo $FILEBASE
	seqret -sequence $FILEBASE.guidance.edit.fasta -outseq $FILEBASE.seqret.phylip -osformat phylip  
	head -q -n1 $FILEBASE.seqret.phylip > $FILEBASE.phy
	tail -q -n +1 $FILEBASE.guidance.edit.fasta >> $FILEBASE.phy   	
	done

