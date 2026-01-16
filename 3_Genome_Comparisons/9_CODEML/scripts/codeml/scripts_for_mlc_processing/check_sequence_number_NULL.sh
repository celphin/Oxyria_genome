#!/bin/bash
FALIST=`ls -1 *.mlc`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=$(echo $FAFILE | cut -f1 -d.)
        echo $FILEBASE
	sed -n '3{p;q}' $FILEBASE.fixed.mlc > $FILEBASE.line3
	awk '{print $1}' $FILEBASE.line3 > $FILEBASE.NoSeq.txt
        done
