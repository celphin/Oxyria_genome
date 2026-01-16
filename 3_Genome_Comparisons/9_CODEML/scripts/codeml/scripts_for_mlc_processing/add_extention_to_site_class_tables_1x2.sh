#!/bin/bash
FALIST=`ls -1 *.site.classes.txt`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=$(echo $FAFILE | cut -f1 -d.)
        echo $FILEBASE
	mv $FILEBASE.site.classes.txt $FILEBASE.1x2.site.classes.txt
        done
