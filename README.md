


several steps to the process

1. unzip and re-sort the source illunima csv from sample-chrom-position to chrom-position-sample
2. convert the sorted file into vcf
3. polish the vcf and index

this is necessary because the vcf has multiple samples on each row


prepare input

get the header that won't be sorted:

```sh
unzip -p FinalReport.zip | head | gzip > input.csv.gz
```

get and sort the rest, removing reference sample on the way.
- this should be sorted correctly, except that chrMT will come before chrX, chrY and chrXY
- not sure if that's going to cause tabix to break or not

```sh
unzip -p FinalReport.zip | tail -n+11 | grep -v NA12878 |sort -Vt , -k 9 -k 10 -k 4 -k 5 | gzip >> input.csv.gz
```

Running the script is done like:

```sh
gunzip -c input.txt.gz | python illumina2vcf --fasta Homo_sapiens_assembly38.fasta --tab | bgzip > output.vcf.gz
```

- these steps can be piped together
```sh
cat <(unzip -p FinalReport.zip | head; \
unzip -p FinalReport.zip | tail -n+11 | \
grep -v NA12878 |sort -Vt , -k 9 -k 10 -k 4 -k 5) | \
python illumina2vcf --fasta Homo_sapiens_assembly38.fasta | \
bgzip > out.vcf.gz
```

Polishing includes resorting (chrX SNPs out of order because chrXY get converted to chrX) and tabix indexing

```sh
bcftools sort out.vcf.gz | bcftools view -s ^NA12878 --no-version --no-update -Oz > out.clean.vcf.gz
tabix out.clean.vcf.gz
```
blocklist
---------

SNPs listed in a file can be excluded using the `--blocklist` option

GSAMD-24v3-0-EA_20034606_A2_blacklist.txt contains all SNPs on the GSA+MD chip that have
a call rate <95% across 267 LRRK2 samples processed between mid Jan and mid June 2022
plus 1475 GOLD samples. Plus SNPs that have HWE violations < e-6 in either dataset
(except for chrX where it was only calculated for female samples from the GOLD study).
There are a number of SNPs that have large allele frequency deviations relative to 1KGP
but have not (yet) included them in this list.

development
===========

```sh
python3 -m venv venv # Create virtual environment
source venv/bin/activate # Activate virtual environment
pip install -e '.[dev]'  # Install using pip including development extras
pre-commit install  # Enable pre-commit hooks
pre-commit run --all-files  # Run pre-commit hooks without committing
# Note pre-commit is configured to use:
# - seed-isort-config to better categorise third party imports
# - isort to sort imports
# - black to format code
pip-compile  # Freeze dependencies
pytest  # Run tests
coverage run -m pytest && coverage report -m  # Run tests, print coverage
mypy .  # Type checking
pipdeptree  # Print dependencies
```

Global git ignores per https://help.github.com/en/github/using-git/ignoring-files#configuring-ignored-files-for-all-repositories-on-your-computer

For release to PyPI see https://packaging.python.org/tutorials/packaging-projects/

To add a new development dependency, add to `requirements-dev.in` then run `pip-compile requirements.in`

To add a new dependency, add to `requirements.in` then run `pip-compile`
