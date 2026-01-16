#!/bin/bash
# Job name: grep_commands
FALIST=`ls -1 *.dist.sorted`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=${FAFILE%.dist.sorted}
        echo $FILEBASE
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Ficmic") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Cansat") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Distri") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Malbac") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Pyrbet") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt      
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Erijap") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Prunavi") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Prunper") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk  '{ printf( "%s ", $1 ); } END { printf( "\n" ); }' $FILEBASE.oneofeach.txt > $FILEBASE.oneofeach.transposed.txt
	done
