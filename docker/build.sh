#!/bin/bash
set -e
set -x

docker build \
  -t illumina2vcf:latest \
  -f docker/Dockerfile \
  . 