#!/bin/bash
set -e
set -x

# build the slim default version
docker build \
  -t illumina2vcf:latest \
  -f docker/Dockerfile \
  .

# build the reference version (~3GB bigger)
docker build \
  -t illumina2vcf:latest-ref \
  -f docker/Dockerfile.ref \
  .
