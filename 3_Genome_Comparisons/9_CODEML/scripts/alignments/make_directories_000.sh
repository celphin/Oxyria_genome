#!/bin/sh
i=0; 
for f in *; 
do 
    d=dir.$(printf %03d $((i/100+1))); 
    mkdir -p $d; 
    mv "$f" $d; 
    let i++; 
done
