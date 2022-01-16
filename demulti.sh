#!/bin/bash

#SBATCH --job-name=barcode
#SBATCH --account="" 
#SBATCH --cpus-per-task=48 
#SBATCH --nodes=1
#SBATCH --time=12:00:00 
#SBATCH --partition="ct560" # or ctest if time<0.5 hr
#SBATCH --mem=96G
#SBATCH -o %j.log           #path to std output
#SBATCH -e %j.err           #path to std err 

TEMP=/work/ylllab2021/temp

# mkdir $TEMP/LIB/fastq/split_pool; mkdir $TEMP/LIB/fastq/split_pool/pass

# rm -rf $TEMP/LIB/fastq/split_1/fail 
# rm -rf $TEMP/LIB/fastq/split_2/fail 
# rm -rf $TEMP/LIB/fastq/split_3/fail 
# rm -rf $TEMP/LIB/fastq/split_4/fail 

# cp $TEMP/LIB/fastq/split_1/pass/*.fastq $TEMP/LIB/fastq/split_pool/pass/
# cp $TEMP/LIB/fastq/split_2/pass/*.fastq $TEMP/LIB/fastq/split_pool/pass/
# cp $TEMP/LIB/fastq/split_3/pass/*.fastq $TEMP/LIB/fastq/split_pool/pass/
# cp $TEMP/LIB/fastq/split_4/pass/*.fastq $TEMP/LIB/fastq/split_pool/pass/

mkdir $TEMP/LIB/demulti

~/app/ont-guppy-cpu_v5/bin/guppy_barcoder -r --require_barcodes_both_ends -i $TEMP/LIB/fastq/ -s $TEMP/LIB/demulti --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg" 


cp $TEMP/LIB/fastq/split_1/sequencing_summary.txt $TEMP/LIB/fastq/

sed '1d' $TEMP/LIB/fastq/split_2/sequencing_summary.txt >> $TEMP/LIB/fastq/sequencing_summary.txt
sed '1d' $TEMP/LIB/fastq/split_3/sequencing_summary.txt >> $TEMP/LIB/fastq/sequencing_summary.txt
sed '1d' $TEMP/LIB/fastq/split_4/sequencing_summary.txt >> $TEMP/LIB/fastq/sequencing_summary.txt
