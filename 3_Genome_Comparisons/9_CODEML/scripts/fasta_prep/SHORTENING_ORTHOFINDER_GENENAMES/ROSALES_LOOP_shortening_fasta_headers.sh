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
	sed -i 's/Cansat_SHORT_pep_//g' $FILEBASE.fa #species1
	sed -i 's/Distri_SHORT_pep_//g' $FILEBASE.fa #species2
	sed -i 's/Erijap_SHORT_pep_//g' $FILEBASE.fa #species3
	sed -i 's/Ficmic_SHORT_pep_//g' $FILEBASE.fa #species4
	sed -i 's/Malbac_SHORT_pep_//g' $FILEBASE.fa #species5
	sed -i 's/Prunper_SHORT_pep_//g' $FILEBASE.fa #species6
	sed -i 's/Pyrbet_SHORT_pep_//g' $FILEBASE.fa #species7
	sed -i 's/Prunavi_SHORT_pep_//g' $FILEBASE.fa #species8
	done
