#!/bin/bash
# Job name: Rearrange distmat output files
FALIST=`ls -1 *.distmat`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=${FAFILE%.distmat}
        echo $FILEBASE
        sed 1,7d $FILEBASE.distmat > $FILEBASE.distmat.noheader #removing the seven first lines
        awk '{for (i=NF;i>0;i--){printf $i" "};printf "\n"}' $FILEBASE.distmat.noheader > $FILEBASE.reversed
	awk '{print $2,$3}' $FILEBASE.reversed > $FILEBASE.two
        sed '/nan/d' $FILEBASE.two > $FILEBASE.noNAs
	sort -n -k 2 $FILEBASE.noNAs > $FILEBASE.dist.sorted 
	done
