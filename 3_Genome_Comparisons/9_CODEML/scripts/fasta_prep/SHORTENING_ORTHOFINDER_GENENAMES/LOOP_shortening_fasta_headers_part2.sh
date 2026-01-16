#!/bin/bash
#Description: Script for shortening fasta headers from orthofinder
#This will delete everything after a pipe
######################################################################################
FALIST=`ls -1 *.fa`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=$(echo $FAFILE | cut -f1 -d.)
        echo $FILEBASE
	sed -i 's@|[^ ]*@@g' $FILEBASE.fa
	done
