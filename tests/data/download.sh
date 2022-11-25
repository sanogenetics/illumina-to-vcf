#!/bin/bash
set -x
set -e
#set -o pipefail

# Homo_sapiens_assembly38.trim.fasta
samtools faidx "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta" 'chr1:1-3000000' 'chr2:1-3000000' 'chr10:1-3000000' 'chrX:1-3000000' 'chrY:1-3000000' |
    sed 's/:1-[0-9]*//' >tests/data/Homo_sapiens_assembly38.trim.fasta
# Homo_sapiens_assembly38.trim.fasta.fai
samtools faidx tests/data/Homo_sapiens_assembly38.trim.fasta

aws s3 cp s3://sano-public/GSA-24v3-0_A2.csv.gz tests/data/GSA-24v3-0_A2.csv.gz
gunzip -c tests/data/GSA-24v3-0_A2.csv.gz | head -n 8 >tests/data/GSA-24v3-0_A2.trim.csv
gunzip -c tests/data/GSA-24v3-0_A2.csv.gz |
    grep -E '(,38,1,)|(,38,2,)|(,38,10,)|(,38,X,)|(,38,Y,)' |
    awk -F',' '$11<3000000' >>tests/data/GSA-24v3-0_A2.trim.csv
gzip -c tests/data/GSA-24v3-0_A2.trim.csv >tests/data/GSA-24v3-0_A2.trim.csv.gz
rm tests/data/GSA-24v3-0_A2.csv.gz
