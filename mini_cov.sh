#!/bin/bash
set -e

# data
run_name="MMDDMN00X-nnBC"
barcode=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24")
threads=24

primerSchemes="/home/ylllab2021/app/artic-ncov2019/primer_schemes"
fastq_dir="/work/ylllab2021/temp/LIB/demulti"
fast5_dir="/work/ylllab2021/temp/LIB/fast5_pass"
seq_summary_dir="/work/ylllab2021/temp/LIB/fastq/sequencing_summary.txt"

# env

source activate artic 

# pipeline commands

for i in ${barcode[@]}

do
cm_guppyplex="artic guppyplex \
              --min-length 100 \
              --max-length 30000 \
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

  samtools mpileup -d 100000 -uf ${primerSchemes}/nCoV-2019/V3/nCoV-2019.reference.fasta ${run_name}${i}.sorted.bam | bcftools call -cv -Oz -o ${run_name}${i}.pile.gz

  tabix ${run_name}${i}.pile.gz

  bcftools consensus -f ${primerSchemes}/nCoV-2019/V3/nCoV-2019.reference.fasta ${run_name}${i}.pile.gz -i '(type="snp")&((DP4[0]+DP4[1])<(DP4[2]+DP4[3]))' > ${run_name}${i}.raw.fasta
 
  sed -i 's/>.*/>pileup/' ${run_name}${i}.raw.fasta
  
  cat ${run_name}${i}.muscle.in.fasta ${run_name}${i}.raw.fasta > ${run_name}${i}.com.fasta

  muscle -in ${run_name}${i}.com.fasta -out ${run_name}${i}.align.fasta

  mkdir ./temp 

  mv *'barcode'${i}* ./temp
  mv ${run_name}${i}* ./temp

  cd ./temp
  multiqc . 
  cd ..

  mv ./temp ${run_name}${i}

  echo "######################################    barcode${i} done"

done

