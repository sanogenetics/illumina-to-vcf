import argparse
import logging
import sys

from illumina2vcf import Converter

logger = logging.getLogger(__name__)

if __name__ == "__main__":

    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
    parser = argparse.ArgumentParser()
    # TODO make non-optional argument
    parser.add_argument("--fasta", help="path to reference fasta, and .fai index")
    parser.add_argument(
        "--tab", "-t", action="store_const", const="\t", dest="delim", help="tab as delimitor in source"
    )
    parser.add_argument(
        "--comma", "-c", action="store_const", const=",", dest="delim", help="comma as delimitor in source"
    )
    parser.add_argument(
        "--blocklist", default="ref/GSAMD-24v3-0-EA_20034606_A2_blocklist.txt", help="path to probe name blocklist"
    )
    args = parser.parse_args()

    converter = Converter(args.fasta, args.fasta + ".fai", args.blocklist, args.delim)

    # read from stdin as uncompressed text
    # write to stdout as uncompressed vcf
    # TODO handle zip or gzip compression in
    # TODO handle blockgzip compression out
    converter.convert(sys.stdin, sys.stdout)
