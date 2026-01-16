#!/bin/bash
# Job name: grep_commands
FALIST=`ls -1 *.dist.sorted`
i=1
for FAFILE in $FALIST
	do
	FILEBASE=${FAFILE%.dist.sorted}
	echo $FILEBASE
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Avecar") {print $i;exit;}}}' $FILEBASE.dist.sorted > $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Jatcur") {print $i;exit;}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Kanobo") {print $i;exit;}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Linusi") {print $i;exit;}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Merann") {print $i;exit;}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Poptri") {print $i;exit;}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Salpur") {print $i;exit;}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Poptre") {print $i;exit;}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{ printf( "%s ", $1 ); } END { printf( "\n" ); }' $FILEBASE.oneofeach.txt > $FILEBASE.oneofeach.transposed.txt
	done
