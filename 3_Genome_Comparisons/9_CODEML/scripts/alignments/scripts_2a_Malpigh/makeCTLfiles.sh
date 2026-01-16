#!/bin/bash
#Script to preparing control files for paml runs  
# Usage: makeCTLfiles.sh
# execute in directory with phylip files

FALIST=`ls -1 *.phy`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=$(echo $FAFILE | cut -f1 -d.)
        echo $FILEBASE
        cp BULKRUN.Malpigh.alt.ctl $FILEBASE.alt.ctl
        sed -i -e "s/REPLACE/$FILEBASE/g" $FILEBASE.alt.ctl
        done
