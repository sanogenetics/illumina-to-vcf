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
    parser.add_argument("--manifest", default="", help="path to Bead Pool Manifest file (csv format)")
    args = parser.parse_args()

    if args.fasta.startswith("s3://"):
        converter = Converter(
            fsspec.open(args.fasta), fsspec.open(args.fasta + ".fai"), args.blocklist, args.manifest, args.delim
        )
    else:
        converter = Converter(args.fasta, args.fasta + ".fai", args.blocklist, args.manifest, args.delim)

    # read from stdin as uncompressed text
    # write to stdout as uncompressed vcf
    converter.convert(sys.stdin, sys.stdout)
