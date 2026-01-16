#!/bin/bash
# Job name: grep_commands
FALIST=`ls -1 *.dist.sorted`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=${FAFILE%.dist.sorted}
        echo $FILEBASE
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "DnA") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "DnN") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Dhis") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Dnem") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Acan") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Avern") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Athal") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Aalp") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk  '{ printf( "%s ", $1 ); } END { printf( "\n" ); }' $FILEBASE.oneofeach.txt > $FILEBASE.oneofeach.transposed.txt
	done
