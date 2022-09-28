#!/bin/bash
set -e
set -x
set -u

# arg 1 is the input s3 url
# arg 2 is the output s3 url
# args other are passed to illumina2vcf

# AWS_CLI_ARGS environment variable can be set to e.g. --endpoint-url=http://moto:5000
# tabix has its own variables e.g. HTS_S3_HOST=moto:5000 and HTS_S3_ADDRESS_STYLE=path
# may also need to set AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY even if IAM not implemented
# may also need to enforce HTTP not HTTPS in scheme e.g. s3+http://....

# first have to sort the input rows
# by sample, chromosome, then position _numerically_ i.e. 1,2,10
# do conversion, then sort again
# generate tabix index of output

# check compression via filename
if [[ $1 = *.zip ]]
then
  # use funzip to stream unzip without needing the file list at end of zip archive
  DECOMPRESSOR=funzip
elif [[ $1 = *.gz ]]
then
  DECOMPRESSOR=gunzip
else
  echo "Unrecognized file ending $1"
  exit 1
fi

# sort by chromosome, position, probe num, samplenum
{ aws $AWS_CLI_ARGS s3 cp $1 - | $DECOMPRESSOR | head ; aws $AWS_CLI_ARGS s3 cp $1 - | $DECOMPRESSOR | tail -n+11 | sort -V -k9,10 -k4,5 ; } | python3.9 -m illumina2vcf ${@:3} | bcftools sort -Oz | aws $AWS_CLI_ARGS s3 cp - $2

s3role tabix $2
