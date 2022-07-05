# use base image that already has AWS CLI & HTSLIB+SAMTools+VCFTools
FROM sanogenetics/htslib-bcftools-samtools:1.15.1

WORKDIR /app
# install reference genome - ~3GB
# Docker ADD always downloads to verify cache
RUN mkdir ref
RUN curl "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta" \
  -o "ref/Homo_sapiens_assembly38.fasta"
RUN curl "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai" \
  -o "ref/Homo_sapiens_assembly38.fasta.fai"

# install reference extras
COPY ref/GSAMD-24v3-0-EA_20034606_A2_blocklist.txt ref/

# install Python & other dependencies
RUN apt-get update \
  && apt-get install -y \
    python3.9 \
    python3.9-distutils \
    libcurl4-gnutls-dev \
    unzip
ADD https://bootstrap.pypa.io/get-pip.py get-pip.py
RUN python3.9 get-pip.py

# install requirements ahead of main code
# not strictly necessary but uses Docker layer caching for faster repeat building
COPY requirements.txt .
RUN pip install -r requirements.txt

# install the actual code
COPY illumina2vcf ./illumina2vcf
COPY setup.py .
COPY README.md .
COPY requirements-dev.txt .
RUN pip install .

# add the wrapper script
COPY docker/run.sh .

CMD ["sh", "run.sh"]
