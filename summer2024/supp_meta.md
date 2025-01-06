#### 1 remove human

```minimap2 -ax map-ont ~/work/data/references/Homo_sapiens.GRCh38.cdna.all.fa.gz MDA_sample02.fq | samtools fastq -n -f 4 - > MDA_sample02_nonhuman.fq```


#### 2 classify reads 

```kraken2 --db ~/work/data/references/k2_viral_20240904 --thread 4 --minimum-hit-groups 3 --report-minimizer-data --report MDA_sample02.k2report.txt MDA_sample02_nonhuman.fq  > MDA_sample02.kraken.txt```


#### 3 bracken (optional)

```bracken -d ~/work/data/references/k2_viral_20240904 -r 200 -i MDA_sample02.k2report.txt -l S -o MDA_sample02.bracken.tsv```


#### 4 extract species reads 

extract_kraken_reads.py is available [here](https://github.com/JenniferLu717/KrakenTools)

```python ~/work/KrakenTools/extract_kraken_reads.py -k MDA_sample02.kraken.txt --include-children -s MDA_sample02_nonhuman.fq -t 2955291 -r MDA_sample02.k2report.txt -o MDA_sample02_FLUA.fa```


#### 5 viz mapping 

```
samtools faidx ~/work/data/references/H3_con-genome.fa

minimap2 -ax map-ont ~/work/data/references/H3_con-genome.fa MDA_sample02_FLUA.fa > MDA_sample02_FLUA.sam

samtools view -bS MDA_sample02_FLUA.sam | samtools sort -o MDA_sample02_FLUA.sort.bam - && samtools index MDA_sample02_FLUA.sort.bam
```


#### virus genomes can be prepared by:

```datasets download virus genome taxon 335341 --filename virus.zip```

```
unzip virus.zip

cp ncbi_dataset/data/genomic.fna .

grep -v "^>" genomic.fna > sequences_only.txt

tr -d '\n' < sequences_only.txt > concatenated_sequence.txt

echo ">taxid_335341_genome" | cat - concatenated_sequence.txt > taxid_335341_genome.fa
```

note: 
1. download genome sequences at subspecies levels (```335341``` = H3N2 subtype) unless you known how many avaialble seuqneces in the species levels 

2. ```datasets``` is the function of _ncbi-datasets-cli_. Install by ```conda install -c conda-forge ncbi-datasets-cli```.











