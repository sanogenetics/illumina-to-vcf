


several steps to the process

1. unzip and re-sort the source illunima csv from sample-chrom-position to chrom-position-sample
2. convert the sorted file into vcf
3. polish the vcf and index

this is necessary because the vcf has multiple samples on each row


prepare input

get the header that won't be sorted:

```sh
unzip -p ~/sano/playground/eurofins_ftp_us/2022_04_19/N23_1_FinalReport.zip | head | gzip > input.header.csv.gz
```

get and sort the rest, removing reference sample on the way:

```sh
unzip -p ~/sano/playground/eurofins_ftp_us/2022_04_19/N23_1_FinalReport.zip | tail -n+11 | grep -v NA12878 |sort -t , -k 9 -k 10 -k 4 -k 5 | gzip > input.body.csv.gz
```

now combine them:

```sh
gunzip -c input.header.csv.gz input.body.csv.gz | gzip > input.csv.gz
```

Running the script is done like:

```sh
gunzip -c input.csv.gz | python run.py | bgzip > out.vcf.gz
```

Polishing includes resorting (sort command earlier uses text sort, so position 1000 is before 10) and tabix indexing

```sh
bcftools sort out.vcf.gz | bcftools view -s ^NA12878 --no-version --no-update -Oz > out.clean.vcf.gz
tabix out.clean.vcf.gz
```