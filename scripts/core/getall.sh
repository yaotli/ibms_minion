mkdir coverage
mkdir bam
mkdir species

cp result/*/*coverage.txt coverage
cp result/*/*.bam bam
 ls -l result | awk '{print $9}' | xargs -iR cp result/R/k2.report_bracken_species.txt species/R.txt
