Illumina2VCF
===========

This is tool for converting Illumina FinalReport microarray results into VCF files. This reconciles multiple probes for the same location, as well as including the reference allele appropriately.

usage
-----

several steps to the process

1. unzip and re-sort the source illunima csv from sample-chrom-position to chrom-position-sample
2. convert the sorted file into vcf
3. polish the vcf and index

this is necessary because the vcf has multiple samples on each row

get the header that won't be sorted:

```sh
unzip -p FinalReport.zip | head | gzip > input.tsv.gz
```

get and sort the rest, removing reference sample on the way.
- this should be sorted correctly, except that chrMT will come before chrX, chrY and chrXY

```sh
unzip -p FinalReport.zip | tail -n+11 | grep -v NA12878 | sort -Vt , -k 9 -k 10 -k 4 -k 5 | gzip >> input.tsv.gz
```

Running the script is done like:

```sh

gunzip -c input.tsv.gz | python3 -m illumina2vcf --manifest GSA-24v3-0_A2.csv --tab Homo_sapiens_assembly38.fasta | bgzip > output.vcf.gz

```

Polishing includes resorting (chrX SNPs out of order because chrXY get converted to chrX) and tabix indexing

```sh
bcftools sort output.vcf.gz | bcftools view -s ^NA12878 --no-version --no-update -Oz > out.clean.vcf.gz
tabix out.clean.vcf.gz
```

Note these steps can be piped together
```sh
{ unzip -p FinalReport.zip | head; \
  unzip -p FinalReport.zip | tail -n+11 \
  | grep -v NA12878 | sort -Vt , -k 9 -k 10 -k 4 -k 5 ; \
 } | python3 -m illumina2vcf --manifest GSA-24v3-0_A2.csv --tab Homo_sapiens_assembly38.fasta \
| bcftools sort -Oz > out.vcf.gz
```

manifest
---------
[Bead Pool Manifest csv file](https://emea.support.illumina.com/bulletins/2016/05/infinium-genotyping-manifest-column-headings.html)
with information about the sequence of indels (optional; indels will be excluded from vcf if not provided)

blocklist
---------

SNPs listed in a file can be excluded using the `--blocklist` option

GSAMD-24v3-0-EA_20034606_A2_blacklist.txt contains all SNPs on the GSA+MD chip that have
a call rate <95% across 267 LRRK2 samples processed between mid Jan and mid June 2022
plus 1475 GOLD samples. Plus SNPs that have HWE violations < e-6 in either dataset
(except for chrX where it was only calculated for female samples from the GOLD study).
There are a number of SNPs that have large allele frequency deviations relative to 1KGP
but have not (yet) included them in this list.

reference
---------

To correctly identify the reference allele, illumina2vcf needs to access a reference genome. This must be in
uncompressed fasta format and have an accompanying .fai index file too.

A build 38 reference file (3GB in size) is available from https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta and has been copied to s3://sano-public/Homo_sapiens_assembly38.fasta in the `latest-ref` tagged container.

For testing, the full reference is too large (~3GB). So there is an alternative that uses a subset of the reference tagged `latest-test-ref` that is ~5MB and limited to:
 - chr1,2,10 1MB for chromosome ordering and autosomal testing
 - chrX,Y 3MB for sex chromosome and PAR testing
 - chrM for mitochondial chromosome testing


S3 access is supported via the [s3fs](http://s3fs.readthedocs.io/en/latest/) library.

Note: uses PyFaidx to read the fasta reference. In theory, this supports a compressed reference. However, to
do so it decompressess all of a chromosome at a time, and my default does this for every access without
caching. So for the purposes of illumina2vcf it is much faster to use a decompressed copy directly.

development
===========

```sh
# Create virtual environment
python3 -m venv venv
# Activate virtual environment
source venv/bin/activate

# first install pip tools
pip install pip-tools
# Install exact requirements
pip-sync requirements.txt requirements-dev.txt
# Install editable with development extras
pip install -e '.[dev]'

# before making any commits
# Enable pre-commit hooks
pre-commit install
# Run pre-commit hooks without committing
pre-commit run --all-files
# Note pre-commit is configured to use:
# - isort to sort imports
# - black to format code

# before making any code changes
# Run tests
pytest
# Run tests, print coverage
coverage run -m pytest && coverage report -m
# Type checking
mypy illumina2vcf

# after making any dependency changes
# Freeze dependencies
pip-compile
# Freeze dev dependencies
pip-compile requirements-dev.in
# then run pip-sync again to install them!
# Print dependencies
pipdeptree

```

Global git ignores per https://help.github.com/en/github/using-git/ignoring-files#configuring-ignored-files-for-all-repositories-on-your-computer

For release to PyPI see https://packaging.python.org/tutorials/packaging-projects/

To add a new development dependency, add to `requirements-dev.in` then run `pip-compile requirements-dev.in`

To add a new dependency, add to `requirements.in` then run `pip-compile`


docker
------

To build this container, use a command like:

```
docker/build.sh
```

Note: `--rm` means to remove intermediate containers after a build. You may want to omit this if developing locally to utilize docker layer caching.

Note: `--progress=plain` may be useful to see all intermediate step logs.

Push to AWS ECR with:

```
aws ecr get-login-password --region eu-west-2 | docker login --username AWS --password-stdin 244834673510.dkr.ecr.eu-west-2.amazonaws.com
docker tag illumina2vcf:latest 244834673510.dkr.ecr.eu-west-2.amazonaws.com/illumina2vcf:latest
docker push 244834673510.dkr.ecr.eu-west-2.amazonaws.com/illumina2vcf:latest
docker tag illumina2vcf:latest-ref 244834673510.dkr.ecr.eu-west-2.amazonaws.com/illumina2vcf:latest-ref
docker push 244834673510.dkr.ecr.eu-west-2.amazonaws.com/illumina2vcf:latest-ref
docker tag illumina2vcf:latest-test-ref 244834673510.dkr.ecr.eu-west-2.amazonaws.com/illumina2vcf:latest-test-ref
docker push 244834673510.dkr.ecr.eu-west-2.amazonaws.com/illumina2vcf:latest-test-ref
docker logout
```

Push to DockerHub with:

```
docker login --username sanogenetics
docker tag illumina2vcf:latest sanogenetics/illumina2vcf:latest
docker push sanogenetics/illumina2vcf:latest
docker tag illumina2vcf:latest-ref sanogenetics/illumina2vcf:latest-ref
docker push sanogenetics/illumina2vcf:latest-ref
docker tag illumina2vcf:latest-test-ref sanogenetics/illumina2vcf:latest-test-ref
docker push sanogenetics/illumina2vcf:latest-test-ref
docker logout
```

Note: will prompt for password.
