import argparse
import logging
import sys

import fsspec

from . import Converter

logger = logging.getLogger(__name__)

if __name__ == "__main__":

    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="path to reference fasta, and .fai index")
    parser.add_argument(
        "--tab", "-t", action="store_const", const="\t", dest="delim", help="tab as delimitor in source"
    )
    parser.add_argument(
        "--comma", "-c", action="store_const", const=",", dest="delim", help="comma as delimitor in source"
    )
    parser.add_argument("--blocklist", default="", help="path to probe name blocklist")
    args = parser.parse_args()

    reference = args.fasta
    reference_fai = reference + ".fai"
    if reference.startswith("s3://"):
        reference = fsspec.open(reference)
        reference_fai = fsspec.open(reference_fai)

    converter = Converter(reference, reference_fai, args.blocklist, args.delim)

    # read from stdin as uncompressed text
    # write to stdout as uncompressed vcf
    converter.convert(sys.stdin, sys.stdout)
