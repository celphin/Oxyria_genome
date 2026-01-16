#!/bin/sh
FALIST=`ls -1 *.fasta`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=${FAFILE%.fasta}
        echo $FILEBASE
        mafft --auto $FILEBASE.fasta > $FILEBASE.aligned.fasta
        done

