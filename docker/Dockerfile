# first build a wheel file
FROM ubuntu:mantic as buildenv

WORKDIR /app

# install Python
RUN apt-get update \
  && apt-get install -y python3.11 python3-pip

RUN pip install hatch --break-system-packages

# install the actual code
COPY illumina2vcf ./illumina2vcf
COPY pyproject.toml .
COPY README.md .

RUN hatch build -t wheel


# use base image that already has AWS CLI & HTSLIB+SAMTools+VCFTools
FROM htslib-bcftools-samtools:latest

WORKDIR /app
# install Python
RUN apt-get update \
  && apt-get install -y python3.11 python3-pip

# install the wheel file built previously
COPY --from=buildenv /app/dist/illumina2vcf-0.0.1-py3-none-any.whl /app
RUN pip install /app/illumina2vcf-0.0.1-py3-none-any.whl --break-system-packages \
  && rm /app/illumina2vcf-0.0.1-py3-none-any.whl

# add the wrapper script
ENV AWS_CLI_ARGS=
COPY docker/run.sh .

CMD ["sh", "run.sh"]
