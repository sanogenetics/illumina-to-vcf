


several steps to the process

1. unzip and re-sort the source illunima csv from sample-chrom-position to chrom-position-sample
2. convert the sorted file into vcf
3. polish the vcf and index

this is necessary because the vcf has multiple samples on each row


prepare input

get the header that won't be sorted:

```sh
unzip -p ~/sano/playground/eurofins_ftp_us/2022_04_19/N23_1_FinalReport.zip | head | gzip > input.csv.gz
```

get and sort the rest, removing reference sample on the way.
- this should be sorted correctly, except that chrMT will come before chrX, chrY and chrXY
- not sure if that's going to cause tabix to break or not

```sh
unzip -p ~/sano/playground/eurofins_ftp_us/2022_04_19/N23_1_FinalReport.zip | tail -n+11 | grep -v NA12878 |sort -Vt , -k 9 -k 10 -k 4 -k 5 | gzip >> input.csv.gz
```

Running the script is done like:

```sh
gunzip -c input.csv.gz | python run.py | bgzip > out.vcf.gz
```

- these steps can be piped together (also removing indels so they don't get logged as errors)
```
grep -v "\[D/I\]" <(unzip -p ../Neuron23/040522-N23_1_MG_FinalReport2.zip | head; \
  unzip -p ../Neuron23/040522-N23_1_MG_FinalReport2.zip | tail -n+11 | \
  grep -v NA12878 |sort -Vt , -k 9 -k 10 -k 4 -k 5) | \
  python run.py --fasta ~/Documents/Sano/local_analyses/QC/Homo_sapiens_assembly38.fasta | \
  bgzip > data/test.vcf.gz
```

```
Polishing includes resorting (sort command earlier uses text sort, so position 1000 is before 10) and tabix indexing

```sh
bcftools sort out.vcf.gz | bcftools view -s ^NA12878 --no-version --no-update -Oz > out.clean.vcf.gz
tabix out.clean.vcf.gz
```