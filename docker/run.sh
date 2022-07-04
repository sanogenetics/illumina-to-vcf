#!/bin/bash
set -e
set -u

aws s3 cp $1 - | bgzip -d | python3.9 -m illumina2vcf ${@:3} | bgzip | aws s3 cp - $2