#!/bin/bash
FALIST=`ls -1 *.alt.mlc`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=$(echo $FAFILE | cut -f1 -d.)
        echo $FILEBASE
        x=$(grep "lnL" $FILEBASE.alt.mlc)
	a=$(grep "proportion" $FILEBASE.alt.mlc)
        b=$(grep "background w" $FILEBASE.alt.mlc)
        c=$(grep "foreground w" $FILEBASE.alt.mlc)
        echo $x > $FILEBASE.site.classes.txt
	echo $a $b $c >> $FILEBASE.site.classes.txt
        done
