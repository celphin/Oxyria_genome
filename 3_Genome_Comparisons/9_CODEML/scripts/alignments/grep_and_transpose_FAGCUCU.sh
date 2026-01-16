#!/bin/bash
# Job name: grep_commands
FALIST=`ls -1 *.dist.sorted`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=${FAFILE%.dist.sorted}
        echo $FILEBASE
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Momcha") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Cucarg") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Cucsat") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Querob") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Betpen") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt      
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Alnglu") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Carill") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Jugreg") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk  '{ printf( "%s ", $1 ); } END { printf( "\n" ); }' $FILEBASE.oneofeach.txt > $FILEBASE.oneofeach.transposed.txt
	done
