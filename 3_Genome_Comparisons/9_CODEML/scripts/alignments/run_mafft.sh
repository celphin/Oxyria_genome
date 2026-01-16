#!/bin/sh
FALIST=`ls -1 *.fa`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=${FAFILE%.fa}
        echo $FILEBASE
        mafft --auto $FILEBASE.fa > $FILEBASE.aligned.fa
        done

