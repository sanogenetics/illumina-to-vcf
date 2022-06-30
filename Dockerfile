FROM python:3.9-slim-buster

WORKDIR /app

# install reference sequence ~3GB ?
#ADD https://sano-public.s3.eu-west-2.amazonaws.com/Homo_sapiens_assembly38.fasta ref/Homo_sapiens_assembly38.fasta
#ADD https://sano-public.s3.eu-west-2.amazonaws.com/Homo_sapiens_assembly38.fasta.fai ref/Homo_sapiens_assembly38.fasta.fai


# install requirements ahead of main code
# not strictly necessary but uses Docker layer caching for faster repeat building
COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . .
RUN pip3 install .

CMD [ "python3", "-m" , "illumina2vcf", "--fasta", "s3://sano-public/Homo_sapiens_assembly38.fasta"]
