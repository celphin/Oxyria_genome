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
	sed -i 's/Avecar_SHORT_pep_//g' $FILEBASE.fa #species1
	sed -i 's/Jatcur_SHORT_pep_//g' $FILEBASE.fa #species2
	sed -i 's/Kanobo_SHORT_pep_//g' $FILEBASE.fa #species3
	sed -i 's/Linusi_SHORT_pep_//g' $FILEBASE.fa #species4
	sed -i 's/Merann_SHORT_pep_//g' $FILEBASE.fa #species5
	sed -i 's/Poptri_SHORT_pep_//g' $FILEBASE.fa #species6
	sed -i 's/Salpur_SHORT_pep_//g' $FILEBASE.fa #species7
	sed -i 's/Poptre_SHORT_pep_//g' $FILEBASE.fa #species8
	done
