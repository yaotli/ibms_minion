#!/bin/bash
set -e

#arg1 
#folder (e.g. LIB/fast5_pass) 
#arg2
#number (e.g. 25 if number of fast5 = 99; >= N/4)

dir=$1
n=$2

cd $dir

mkdir split_1 split_2 split_3 split_4

mv `ls *.fast5 | head -$n` split_1
mv `ls *.fast5 | head -$n` split_2
mv `ls *.fast5 | head -$n` split_3
mv `ls *.fast5 | head -$n` split_4

cd - 

