#!/bin/bash
set -e

# data
run_name="0628MN002-14BC"
barcode=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14")
threads=1

primerSchemes="/home/users/nus/e0081754/app/artic-ncov2019/primer_schemes"
fastq_dir="./demulti"
fast5_dir="./fast5_pass"
seq_summary_dir="./sequencing_summary_FAP76995_7a9c076c.txt"

# env

source activate artic-ncov2019 

# pipeline commands

for i in ${barcode[@]}

do
cm_guppyplex="artic guppyplex \
              --min-length 400 \
              --max-length 700 \
              --directory ${fastq_dir}/barcode${i} \
              --prefix ${run_name}"

cm_minion="artic minion \
           --normalise 500 \
           --threads ${threads} \
           --scheme-directory ${primerSchemes} \
           --read-file ${run_name}_barcode$i.fastq \
           --fast5-directory ${fast5_dir} \
           --strict \
           --sequencing-summary ${seq_summary_dir}
           nCoV-2019/V3 \
           ${run_name}${i}"

function runCml {
    echo 
    "$@"
    echo
    return 
}

  runCml $cm_guppyplex 

  runCml $cm_minion

  mkdir ./temp 

  mv *'barcode'${i}* ./temp
  mv ${run_name}${i}* ./temp

  cd ./temp
  multiqc . 
  cd ..

  mv ./temp ${run_name}${i}

  echo "######################################    barcode${i} done"

done



