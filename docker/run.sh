#!/bin/bash
set -e
set -u

# first have to sort the input rows 
# by sample, chromosome, then position _numerically_ i.e. 1,2,10
# after conversion, sort again (?)

{ aws s3 cp $1 - | funzip | head; \
  aws s3 cp $1 - | funzip | tail -n+11 \
  | grep -v NA12878 | sort -Vt , -k 9 -k 10 -k 4 -k 5
} | python3.9 -m illumina2vcf ${@:3} | bcftools sort -Oz | aws s3 cp - $2

s3role tabix $2
