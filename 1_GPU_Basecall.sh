#!/bin/bash

set -e

TEMP=/work/ylllab2021/temp

echo 'make sure you...'
echo '1 log in ln01.twcc.ai'
echo '2 make fast5_pass.tar.gz is directly under temp' 
echo '3 no other folder under temp called LIB'
echo ''
echo 'press Ctl+c if NOT' 

sleep 5 

echo 'start unzip... it may take awile in TWCC'
tar -xvf $TEMP/fast5_pass.tar.gz -C $TEMP
mkdir -p $TEMP/LIB
mv $TEMP/fast5_pass $TEMP/LIB

echo 'start basecalling'
sbatch basecall_gpu.sh


