#!/bin/bash

#SBATCH --job-name=split_1
#SBATCH --account="" 
#SBATCH --cpus-per-task=48 
#SBATCH --nodes=8
#SBATCH --time=48:00:00 
#SBATCH --partition="ct560" # or ctest if time<0.5 hr
#SBATCH --mem=96G
#SBATCH -o %j.log           #path to std output
#SBATCH -e %j.err           #path to std err 

TEMP=/work/ylllab2021/temp


~/app/ont-guppy-cpu_v5/bin/guppy_basecaller -i $TEMP/LIB/fast5_pass/split_1 -s $TEMP/LIB/fastq/split_1 -q 0 --cpu_threads_per_caller 8 -c dna_r9.4.1_450bps_hac.cfg --num_callers 4 > basecallreport_1.txt

#~/app/ont-guppy-cpu_v5/bin/guppy_barcoder --require_barcodes_both_ends -i ~/temp/test/pass -s ~/temp/demulti --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
