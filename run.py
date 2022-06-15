import argparse
import csv
import itertools
import logging
import sys
import time
from contextlib import contextmanager
from typing import Dict, Generator, List, Tuple

from puretabix.vcf import LINE_START, VCFAccumulator, VCFLine, get_vcf_fsm
from pyfaidx import Fasta

# Sample ID
# GC Score
# SNP Name
# SNP Index
# Sample Index
# Sample Name
# Sample Group
# SNP Aux
# Chr
# Position
# GT Score
# Cluster Sep
# SNP
# ILMN Strand
# Customer Strand
# Top Genomic Sequence
# Plus/Minus Strand
# Allele1 - Plus
# Allele2 - Plus
# Allele1 - Forward
# Allele2 - Forward
# Theta
# R
# X
# Y
# X Raw
# Y Raw
# B Allele Freq
# Log R Ratio
# CNV Value
# CNV Confidence
SAMPLE_ID = "Sample ID"
SNP_NAME = "SNP Name"
ALLELE1 = "Allele1 - Plus"
ALLELE2 = "Allele2 - Plus"

genome_build = "GRCh38"
strandswap = {"A": "T", "T": "A", "C": "G", "G": "C"}

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
logger = logging.getLogger(__name__)


@contextmanager
def timed(message, minimum=1.0):
    start = time.time()
    yield None
    end = time.time()
    duration = end - start
    if duration > minimum:
        logger.info(f"{message} ({duration}s)")


class ConverterError(Exception):
    pass

class DateError(Exception):
    pass

class Converter:
    row_previous_loc = None
    row_previous = None
    vcf_fsm = get_vcf_fsm()
    vcf_accumulator = VCFAccumulator()
    header_done = False
    sample_set = []

    def __init__(self, fasta, fasta_index):
        self.ref = Fasta(fasta)
        self.chromosome_sizes = self.parse_index(fasta_index)
        self.vcf_fsm = get_vcf_fsm()
        self.vcf_accumulator = VCFAccumulator()
        self.sample_set = []

    def parse_index(self, fasta_index):
        chr_lengths = {}
        with open(fasta_index, 'rt') as index_fh:
            index_reader = csv.reader(index_fh, delimiter='\t')
            for line in index_reader:
                chr_lengths[line[0]] = line[1]
        return(chr_lengths)
    
    def ref_lookup(self, chm: str, pos: int) -> str:
        ref_base = str(self.ref[chm][pos-1])
        return ref_base

    def _generate_vcf_header(
        self, input, buildsizes: Dict[str, str], buildname: str
    ) -> Generator[VCFLine, None, None]:
        (date, source) = self._parse_file_header(input)
        # write header
        yield VCFLine.as_comment_key_string("fileformat", "VCFv4.3")
        yield VCFLine.as_comment_key_string("filedate", date)
        yield VCFLine.as_comment_key_string("source", f'"{source}, Sano Genetics"')
        
        # ##FILTER=<ID=PASS,Description="All filters passed">
        yield VCFLine.as_comment_key_dict(
            "FILTER",
            {
                "ID": "PASS",
                "Description": '"All filters passed"'
            } 
        )          

        # ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        yield VCFLine.as_comment_key_dict(

            "FORMAT",
            {
                "ID": "GT",
                "Number": "1",
                "Type": "String",
                "Description": "Genotype",
            }
        )

        # write contig
        for chrom, length in buildsizes.items():
            # ##contig=<ID=1,length=249250621,assembly=GRCh37>
            yield VCFLine.as_comment_key_dict(
                "contig",
                {"ID": chrom, "length": length, "assembly": buildname},
            )


    def _generate_lines(self, input) -> Generator[Dict[str, str], None, None]:
        for row in csv.DictReader(input):
            yield row

    def _parse_file_header(self, input):
        file_header = [row for row in itertools.islice(input, 0, 9)]
        source=file_header[3].split(',')[-1].split('.')[0].lstrip()

        # need some validation on the dates here because there's a good chance they
        # will switch up the format on us at some point
        # this will currently work for month/day/year and year-month-day
        # (I guess also for month-day-year and year/month/day)
        date=file_header[2].split(',')[-1].lstrip().split(' ')[0]
        date_components=date.replace('/', '-').split('-')
        if len(date_components) != 3:
            raise DateError(f"Cannot parse Processing date {date}")
        if len(date_components[2]) == 4:
            date_components = [date_components[2], date_components[0], date_components[1]]
        elif len(date_components[0]) != 4:
            raise DateError(f"Cannot parse Processing date {date}")
        if int(date_components[1]) > 12:
            raise DateError(f"Cannot parse Processing date {date}")
        date_components[1] = date_components[1].zfill(2)
        date_components[2] = date_components[2].zfill(2)
        date=''.join(date_components)
        return(date, source)

    def _generate_line_blocks(
        self, input
    ) -> Generator[List[Dict[str, str]], None, None]:
        block = []
        for row in self._generate_lines(input):
            # is there an existing block that might end?
            if block and (
                block[0]["Chr"] != row["Chr"] or block[0]["Position"] != row["Position"]
            ):
                yield block
                block = []
            block.append(row)
        yield block

    def _generate_vcf_lines(self, input) -> Generator[VCFLine, None, None]:
        for block in self._generate_line_blocks(input):
            assert len(block) >= 1, block
            try:
                with timed(
                    f"block to vcf",
                ):
                    vcfline = self._line_block_to_vcf_line(block)
            except ConverterError as e:
                logger.error(e)
                continue

            if not self.header_done:
                # write the header
                yield VCFLine.as_comment_raw(
                    "\t".join(
                        (
                        "CHROM",
                        "POS",
                        "ID",
                        "REF",
                       "ALT",
                        "QUAL",
                        "FILTER",
                        "INFO",
                        "FORMAT",
                        "\t".join(self.sample_set),  # customizable name of sample
                        )

                    )
                )
                self.header_done = True

            yield vcfline

    def combine_calls(self, previous_calls, new_calls, previous_probes, new_probes):
        if previous_calls == new_calls or previous_calls == new_calls[::-1] or new_calls == ('-', '-'): 
            return(previous_calls)
        elif new_calls[0] == new_calls[1] and new_calls[0] in previous_calls:
            if not new_probes[0] in previous_calls or not new_probes[1] in previous_calls:
                return(previous_calls) # conflicting call is a homozygote where only one of the probes matches the true genotype
            else:
                raise ConverterError("conflict!") # both probes match but different call
        elif previous_calls == ('-', '-'): 
            return(new_calls)
        elif previous_calls[0] == previous_calls[1] and previous_calls[0] in new_calls:
            if new_probes[0] not in previous_probes or new_probes[1] not in previous_probes:
                return(new_calls) # conflicting call is a het where only one of the probes has been tested previously (and previous call is homozygous)
            else:
                raise ConverterError("conflict!") # both probes match but different call
        else:
            raise ConverterError("conflict!") # conflicting hets or conflicting homozygotes
        
    def _simplify_block(self, block: List[Dict[str, str]]):
        combined_probes = {}
        calls = {}
        conflicts = []
        for row in block:
            sampleid = row[SAMPLE_ID]
            strand = row["Plus/Minus Strand"]

            new_probes = list(row["SNP"][1:-1].split("/"))
            # exclude indel probes
            if new_probes[0] not in "ATCG" or new_probes[1] not in "ATCG":
                raise ConverterError(
                    f"Not a SNV probe {row['Chr']}:{row['Position']}"
                )
            if strand == '-':
                new_probes = [strandswap[probe] for probe in new_probes]
            new_calls = (row[ALLELE1], row[ALLELE2])

            if sampleid not in calls:
                calls[sampleid] = new_calls
                combined_probes[sampleid] = set(new_probes)
            else:
                try:
                    calls[sampleid] = self.combine_calls(calls[sampleid], new_calls, combined_probes[sampleid], new_probes)
                except ConverterError as e:
                    #logger.warning(f"{row['Chr']}:{row['Position']}, {sampleid}: {e}")
                    conflicts.append(sampleid)
                combined_probes[sampleid].update(new_probes)

        probed = set().union(*[combined_probes[sampleid] for sampleid in combined_probes.keys()])
        
        # for all samples with conflicting genotype calls, set the call to no call 
        for sampleid in conflicts:
            calls[sampleid] = ("-", "-")

        return (calls, probed)

    def _line_block_to_vcf_line(self, block: List[Dict[str, str]]) -> VCFLine:
        # if block is silly big skip it (this will fail with >100 samples)
        if len(block) > 100 or (self.sample_set and len(block) > 10 * len(self.sample_set)):
            raise ConverterError(
                f"Oversized block {block[0]['Chr']}:{block[0]['Position']}: {len(block)} rows"
            )

        # if we've not got a list of samples yet, get them from this block unfiltered
        if not self.sample_set:
            for row in block:
                if row[SAMPLE_ID] not in self.sample_set:
                    self.sample_set.append(row[SAMPLE_ID])

        # if there are multiple probes that agree, thats fine
        # need to remove conflicting and duplicate rows
        [calls, probed] = self._simplify_block(block)
        #if not calls:
        #    raise ConverterError(
        #        f"Oversimplified block {block[0]['Chr']}:{block[0]['Position']}"
        #    )
        # can assume each row in block has:
        #  same chm, pos
        #  unique sample names
        #  had consistent calls

        snp_names = tuple(sorted(frozenset([r[SNP_NAME] for r in block])))

        chm = block[0]["Chr"]
        # convert pseudoautosomal (XY) to X
        if chm == "XY" or chm == "chrXY":
            chm = "chrX"
        # convert MT to M
        if chm == "MT" or chm == "chrMT":
            chm = "chrM"
        # force chr prefix
        if not chm.startswith("chr"):
            chm = f"chr{chm}"
        if chm not in self.chromosome_sizes.keys():
            raise ConverterError(
                f"Unexpected chromosome {block[0]['Chr']}:{block[0]['Position']}"
            )
        pos = int(block[0]["Position"])

        ref = self.ref_lookup(chm, pos)

        # hande ref/alt split
        if ref not in probed:
            raise ConverterError(
                f"{';'.join(snp_names)}: Reference ({ref}) not probed ({','.join(probed)}) {block[0]['Chr']}:{block[0]['Position']}"
            )
        probed.remove(ref)
        alt = tuple(sorted(probed))
        # alt may not be used, but is what the microarray could check for

        # convert calls
        converted_calls = {}
        for sampleid in calls:
            allele1 = calls[sampleid][0]
            if allele1 not in "ATCG":
                allele1n = "."
            elif allele1 == ref:
                allele1n = 0
            elif allele1 not in alt:
                raise ConverterError(
                    f"Unexpected forward {block[0]['Chr']}:{block[0]['Position']}"
                )
            else:
                allele1n = 1 + alt.index(allele1)

            allele2 = calls[sampleid][1]
            if allele2 not in "ATCG":
                allele2n = "."
            elif allele2 == ref:
                allele2n = 0
            elif allele2 not in alt:
                raise ConverterError(
                    f"Unexpected forward {block[0]['Chr']}:{block[0]['Position']}"
                )
            else:
                allele2n = 1 + alt.index(allele2)

            assert sampleid not in converted_calls
            converted_calls[sampleid] = f"{allele1n}/{allele2n}"

        # always need to have all samples on all rows
        # even if that samples has been filtered out e.g. conflicting probes
        samples = []
        for sampleid in self.sample_set:
            samples.append({"GT": converted_calls.get(sampleid, "./.")})

        vcfline = VCFLine(
            "", "", "", {}, chm, pos, snp_names, ref, alt, ".", ["PASS"], {}, samples
        )
        #vcfline = self.clean_vcf_line(vcfline)

        return vcfline

    def clean_vcf_line(self, vcfline: VCFLine) -> VCFLine:
        # if no samples are called, discard
        if frozenset((s.get("GT", "./.") for s in vcfline.sample)) == frozenset(
            ("./.",)
        ):
            raise ConverterError(f"No calls made {vcfline.chrom}:{vcfline.pos}")

        return vcfline

    def convert(self, input, outfile):
        for headerline in self._generate_vcf_header(input, self.chromosome_sizes, genome_build):
            outfile.write(str(headerline))
            outfile.write("\n")
        for vcfline in self._generate_vcf_lines(input):
            outfile.write(str(vcfline))
            outfile.write("\n")
            # outfile.flush()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", help="path to reference fasta (with fai index)")
    args = parser.parse_args()

    converter = Converter(args.fasta, args.fasta + ".fai")

    # read from stdin
    # write to stdout
    converter.convert(sys.stdin, sys.stdout)
