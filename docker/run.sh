#!/bin/bash
set -e
set -x
set -u

# arg 1 is the input s3 url
# arg 2 is the output s3 url
# args other are passed to illumina2vcf

# first have to sort the input rows
# by sample, chromosome, then position _numerically_ i.e. 1,2,10
# do conversion, then sort again (not needed?)
# generate tabix index of output

# use funzip to stream unzip without needing the file list at end of zip archive
# replace commas to tab incase of comma separated input
# sort by chromosome, position, probe num, samplenum
{ aws ${AWS_CLI_ARGS} s3 cp $1 - | funzip | head; \
  aws ${AWS_CLI_ARGS} s3 cp $1 - | funzip | tail -n+11 \
  | sed 's/\,/\t/g' \
  | sort -V -k9,10 -k4,5
} | python3.9 -m illumina2vcf ${@:3} | bcftools sort -Oz | aws ${AWS_CLI_ARGS} s3 cp - $2

s3role tabix $2
