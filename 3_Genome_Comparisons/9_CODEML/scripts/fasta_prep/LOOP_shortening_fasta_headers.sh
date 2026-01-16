#!/bin/bash
#Description: Script for shortening fasta headers from orthofinder
######################################################################################
FALIST=`ls -1 *.fa`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=$(echo $FAFILE | cut -f1 -d.)
        echo $FILEBASE
	sed -i ':a;N;/^>/M!s/\n//;ta;P;D' $FILEBASE.fa
	sed -i 's/Ampedge_SHORT_pep_//g' $FILEBASE.fa #species1
	sed -i 's/Cercan_SHORT_pep_//g' $FILEBASE.fa #species2
	sed -i 's/Chafas_SHORT_pep_//g' $FILEBASE.fa #species3
	sed -i 's/Dalodo_SHORT_pep_//g' $FILEBASE.fa #species4
	sed -i 's/Faialb_SHORT_pep_//g' $FILEBASE.fa #species5
	sed -i 's/Glymax_SHORT_pep_//g' $FILEBASE.fa #species6
	sed -i 's/Phavul_SHORT_pep_//g' $FILEBASE.fa #species7
	sed -i 's/Vigung_SHORT_pep_//g' $FILEBASE.fa #species8
	done
