#### remove human

```minimap2 -ax map-ont ~/work/data/references/Homo_sapiens.GRCh38.cdna.all.fa.gz MDA_sample02.fq | samtools fastq -n -f 4 - > MDA_sample02_nonhuman.fq```


#### classify reads 

```kraken2 --db ~/work/data/references/k2_viral_20240904 --thread 4 --minimum-hit-groups 3 --report-minimizer-data --report MDA_sample02.k2report.txt MDA_sample02_nonhuman.fq  > MDA_sample02.kraken.txt```


#### bracken (optional)

```bracken -d ~/work/data/references/k2_viral_20240904 -r 200 -i MDA_sample02.k2report.txt -l S -o MDA_sample02.bracken.tsv```


#### extract species reads 

extract_kraken_reads.py is available [here](https://github.com/JenniferLu717/KrakenTools)

```python ~/work/KrakenTools/extract_kraken_reads.py -k MDA_sample02.kraken.txt --include-children -s MDA_sample02_nonhuman.fq -t 2955291 -r MDA_sample02.k2report.txt -o MDA_sample02_FLUA.fa```


#### viz mapping 

```
samtools faidx ~/work/data/references/H3_con-genome.fa

minimap2 -ax map-ont ~/work/data/references/H3_con-genome.fa MDA_sample02_FLUA.fa > MDA_sample02_FLUA.sam

samtools view -bS MDA_sample02_FLUA.sam | samtools sort -o MDA_sample02_FLUA.sort.bam - && samtools index MDA_sample02_FLUA.sort.bam
```

