#!/bin/bash
#Script to make one sorted mlc file 
# execute in directory with all mlc files

find . -name '*.mlc' -exec grep -H 'lnL' {} + > greppomania.file
sed 's:.*/::' greppomania.file > test_greppomania.file
awk '{split($1,a,/:/);$1=a[1]}1' test_greppomania.file > test2_greppomania.file 
awk '{print $1,$5}' test2_greppomania.file > test3_greppomania.file 
awk '{split($1,a,/\./);$1=a[1]}1' test3_greppomania.file > ready_to_sort.file
sort -k1 ready_to_sort.file > SORTED_lnL.file
