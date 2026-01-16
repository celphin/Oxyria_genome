#!/bin/bash
# Job name: grep_commands
FALIST=`ls -1 *.dist.sorted`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=${FAFILE%.dist.sorted}
        echo $FILEBASE
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Dalodo") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Phavul") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Vigung") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Glymax") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Ampedge") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt      
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Chafas") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Faialb") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Cercan") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk  '{ printf( "%s ", $1 ); } END { printf( "\n" ); }' $FILEBASE.oneofeach.txt > $FILEBASE.oneofeach.transposed.txt
	done
