#!/bin/bash
FALIST=`ls -1 *.txt`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=$(echo $FAFILE | cut -f1 -d.)
        echo $FILEBASE
        grep -l -f $FILEBASE.txt ../$FILEBASE*.site.classes.txt >> list_over_appropriate_site_class_tables.file 
        done
