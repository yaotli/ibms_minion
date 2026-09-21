

echo \> Using ${barcode} as Barcode
echo \> Using ${fastq} as FASTQ File Name

nanoq -i demulti/${fastq} -s -t 5 -vvvvv

echo ===

grep ">" core/lib/${reference}

echo ===

# index the reference
samtools faidx core/lib/${reference}

# mapping using minimap2
minimap2 -ax map-ont core/lib/${reference} demulti/${fastq} > barcode${barcode}_aln.sam

# change .sam to .bam
samtools view -bS barcode${barcode}_aln.sam | samtools sort -o barcode${barcode}_sorted.bam - && samtools index barcode${barcode}_sorted.bam

# create an idxstats file
samtools idxstats barcode${barcode}_sorted.bam > barcode${barcode}_sorted_idxstats.txt

echo ===

awk '{if ($3 > max) {max=$3; ref=$1}} END {print ref}' barcode${barcode}_sorted_idxstats.txt

echo ===

samtools view -c -F2308 barcode${barcode}_sorted.bam

echo ===

# sleep 3

if [[ $selected ]]
then
	echo [WARNING] Using manually selecting mode - $selected
	echo $selected | awk '{print $1 * 2 - 1 "," $1 * 2 "p"}'
	echo ---
	echo $selected | awk '{print $1 * 2 - 1 "," $1 * 2 "p"}' | xargs -iR sed -n R core/lib/${reference} > ev_ref.fa
else
	awk 'BEGIN{max=0; t=0; ans=0;} {t=t+1; if ($3 > max) {max=$3; ref=$1; ans=t;}} END {print ans * 2 - 1 "," ans * 2 "p"}' barcode${barcode}_sorted_idxstats.txt
	awk 'BEGIN{max=0; t=0; ans=0;} {t=t+1; if ($3 > max) {max=$3; ref=$1; ans=t;}} END {print ans * 2 - 1 "," ans * 2 "p"}' barcode${barcode}_sorted_idxstats.txt | xargs -iR sed -n R core/lib/${reference} > ev_ref.fa
fi

# get line 5-6 to anothter file called ev_ref
# awk '{if ($3 > max) {max=$3; ref=$1; t=t + 1}} END {print t * 2 - 1 "," t * 2 "p"}' barcode${barcode}_sorted_idxstats.txt
# awk '{if ($3 > max) {max=$3; ref=$1; t=t + 1}} END {print t * 2 - 1 "," t * 2 "p"}' barcode${barcode}_sorted_idxstats.txt | xargs -iR sed -n R core/lib/${reference} > ev_ref.fa

# index the reference
samtools faidx ev_ref.fa

# mapping using minimap2
minimap2 -ax map-ont ev_ref.fa demulti/${fastq} > barcode${barcode}_${name}.sam

# change .sam to .bam
samtools view -bS barcode${barcode}_${name}.sam | samtools sort -o barcode${barcode}_sorted_${name}.bam - && samtools index barcode${barcode}_sorted_${name}.bam

echo ===

bedtools genomecov -d -split -ibam barcode${barcode}_sorted_${name}.bam > barcode${barcode}_coverage.txt

echo ===
echo Use Tablet to see result
echo ===
# sleep 20

# pile up the reads and use bcftools to call the variant
bcftools mpileup -d 50000 -f ev_ref.fa barcode${barcode}_sorted_${name}.bam | bcftools call -cv -Oz -o barcode${barcode}_pile.gz
#samtools mpileup -d 100000 -uf ev_ref.fa barcode${barcode}_sorted_${name}.bam | bcftools call -cv -Oz -o barcode${barcode}_pile.gz

# the variants are saved in a VCF format
tabix barcode${barcode}_pile.gz

# use the variant information to come up a new sequence by editing the reference sequence
bcftools consensus -f ev_ref.fa barcode${barcode}_pile.gz  -i '(type="snp")&((DP4[0]+DP4[1])<(DP4[2]+DP4[3]))' > barcode${barcode}_consensus.fa

echo ===

# map the original .fastq data against human genome
minimap2 -ax map-ont core/lib/host/Homo_sapiens.GRCh38.cdna.all.fa.gz demulti/${fastq} -o human.filtered.sam

# collect unmapped reads in a new fastq
samtools fastq -n -f 4 human.filtered.sam > human.filtered.fastq
pwd
echo ===

# map the original .fastq data against human genome
minimap2 -ax map-ont core/lib/host/Chlorocebus_sabaeus.ChlSab1.1.cdna.all.fa.gz human.filtered.fastq -o monkey.filtered.sam

# collect unmapped reads in a new fastq
samtools fastq -n -f 4 monkey.filtered.sam > monkey.filtered.fastq

echo ===

# map the original .fastq data against human genome
minimap2 -ax map-ont core/lib/host/Aedes_aegypti_lvpagwg.AaegL5.cdna.all.fa.gz monkey.filtered.fastq -o mosquito.filtered.sam

# collect unmapped reads in a new fastq
samtools fastq -n -f 4 mosquito.filtered.sam > human.filtered.fastq

echo ===

# use kraken to map the .fastq to a dataset containing viral species
# k2_viral_20240605 is a folder with indexed database
kraken2 -db ~/k2 --report k2.report.txt human.filtered.fastq > k2.output.log

# summarize the resulg using Bracken
bracken -d k2_viral_20240605 -r 200 -i k2.report.txt  -l S -o bracken.tsv

echo ===

ktImportTaxonomy -t 5 -m 3 k2.report_bracken_species.txt
