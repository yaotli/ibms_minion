#!/bin/bash
set -e

TEMP=/work/ylllab2021/temp

mkdir $TEMP/LIB/fastq

sbatch basecall_1.sh 
sbatch basecall_2.sh
sbatch basecall_3.sh
sbatch basecall_4.sh
