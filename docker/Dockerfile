# use base image that already has AWS CLI & HTSLIB+SAMTools+VCFTools
FROM sanogenetics/htslib-bcftools-samtools:1.15.1

WORKDIR /app
# install Python & other dependencies
RUN apt-get update \
  && apt-get install -y \
    python3.9 \
    python3.9-distutils \
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
ENV AWS_CLI_ARGS=
COPY docker/run.sh .

CMD ["sh", "run.sh"]
