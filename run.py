import argparse
import csv
import itertools
import logging
import sys
import time
from contextlib import contextmanager
from typing import Dict, Generator, List, Tuple

from puretabix.tabix import TabixIndexedFile
from puretabix.vcf import LINE_START, VCFAccumulator, VCFLine, get_vcf_fsm

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

dbsnp_38_chrs = {
    "NC_000001.11": "chr1",
    "NC_000002.12": "chr2",
    "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",
    "NC_000005.10": "chr5",
    "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",
    "NC_000008.11": "chr8",
    "NC_000009.12": "chr9",
    "NC_000010.11": "chr10",
    "NC_000011.10": "chr11",
    "NC_000012.12": "chr12",
    "NC_000013.11": "chr13",
    "NC_000014.9": "chr14",
    "NC_000015.10": "chr15",
    "NC_000016.10": "chr16",
    "NC_000017.11": "chr17",
    "NC_000018.10": "chr18",
    "NC_000019.10": "chr19",
    "NC_000020.11": "chr20",
    "NC_000021.9": "chr21",
    "NC_000022.11": "chr22",
    "NC_000023.11": "chrX",
    "NC_000024.10": "chrY",
}
dbsnp_38_chrs_inv = dict(((j, i) for i, j in dbsnp_38_chrs.items()))
build_38_sizes = {
    "chr1": "248956422",
    "chr2": "242193529",
    "chr3": "198295559",
    "chr4": "190214555",
    "chr5": "181538259",
    "chr6": "170805979",
    "chr7": "159345973",
    "chr8": "145138636",
    "chr9": "138394717",
    "chr10": "133797422",
    "chr11": "135086622",
    "chr12": "133275309",
    "chr13": "114364328",
    "chr14": "107043718",
    "chr15": "101991189",
    "chr16": "90338345",
    "chr17": "83257441",
    "chr18": "80373285",
    "chr19": "58617616",
    "chr20": "64444167",
    "chr21": "46709983",
    "chr22": "50818468",
    "chrX": "156040895",
    "chrY": "57227415",
}

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


class Converter:
    row_previous_loc = None
    row_previous = None
    vcf_fsm = get_vcf_fsm()
    vcf_accumulator = VCFAccumulator()
    header_done = False
    sample_set = []

    def __init__(self, dbsnp_vcf, dbsnp_vcf_index):
        with timed(f"dbsnp index loading"):
            self.dbsnp = TabixIndexedFile.from_files(
                open(dbsnp_vcf, "rb"), open(dbsnp_vcf_index, "rb")
            )
        self.vcf_fsm = get_vcf_fsm()
        self.vcf_accumulator = VCFAccumulator()
        self.sample_set = []

    def ref_alt_lookup(self, chm: str, pos: int) -> Tuple[str, List[str]]:

        # turn the chrX format into other ids

        with timed(f"dbsnp lookup for {chm}:{pos}", 1.0):
            lines = tuple(self.dbsnp.fetch_lines(dbsnp_38_chrs_inv[chm], pos, pos))
        if not lines:
            return "", []
        assert len(lines) == 1, lines

        # sometimes there are multiple lines at the same position
        # this violates VCF spec, but its dbSNPs fault
        # one of these is probably the single nuceotide variant we want
        line = [l + "\n" for l in lines if "VC=SNV" in l][0]

        # parse the line
        self.vcf_fsm.run(line, LINE_START, self.vcf_accumulator)
        vcfline = self.vcf_accumulator.to_vcfline()
        self.vcf_accumulator.reset()

        assert vcfline.pos == pos

        return vcfline.ref, vcfline.alt

    def _generate_vcf_header(
        self, buildsizes: Dict[str, str], buildname: str, samplename: str
    ) -> Generator[VCFLine, None, None]:

        # write header
        yield VCFLine.as_comment_key_string("fileformat", "VCFv4.3")
        # TODO write INFO
        # write SAMPLE
        # write FORMAT
        # ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        yield VCFLine.as_comment_key_dict(
            "FORMAT",
            {
                "ID": "GT",
                "Number": "1",
                "Type": "String",
                "Description": "Genotype",
            },
        )
        # write contig
        for chrom, length in buildsizes.items():
            # ##contig=<ID=1,length=249250621,assembly=GRCh37>
            yield VCFLine.as_comment_key_dict(
                "contig",
                {"ID": chrom, "length": length, "assembly": buildname},
            )

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
                    samplename,  # customizable name of sample
                )
            )
        )

    def _generate_lines(self, input) -> Generator[Dict[str, str], None, None]:
        for row in csv.DictReader(itertools.islice(input, 9, None)):
            yield row

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
            assert len(block) > 1, block
            try:
                with timed(
                    f"block to vcf",
                ):
                    vcfline = self._line_block_to_vcf_line(block)
            except ConverterError as e:
                logger.error(e)
                continue

            if not self.header_done:
                for headerline in self._generate_vcf_header(
                    build_38_sizes, "GRCh38", "\t".join(self.sample_set)
                ):
                    yield headerline
                self.header_done = True

            yield vcfline

    def _simplify_block(self, block: List[Dict[str, str]]) -> List[Dict[str, str]]:

        # drop duplicate lines from a block (same chm,pos,sample,call)
        calls = {}
        samples_to_drop = set()
        for row in block:
            sampleid = row[SAMPLE_ID]
            # TODO handle probes of AA + --
            if sampleid not in calls:
                calls[sampleid] = (row[ALLELE1], row[ALLELE2])
            elif calls[sampleid] != (row[ALLELE1], row[ALLELE2]):
                # previous call was different, drop this sample
                samples_to_drop.add(sampleid)

        rows_keep = {}
        for row in block:
            sampleid = row[SAMPLE_ID]
            if sampleid not in samples_to_drop and sampleid not in rows_keep:
                rows_keep[sampleid] = row

        newblock = list(rows_keep.values())
        return newblock

    def _line_block_to_vcf_line(self, block: List[Dict[str, str]]) -> VCFLine:
        # if block is silly big skip it
        if len(block) > 100:
            raise ConverterError(
                f"Oversized block {block[0]['Chr']}:{block[0]['Position']}"
            )

        # if we've not got a list of samples yet, get them from this block unfiltered
        if not self.sample_set:
            for row in block:
                if row[SAMPLE_ID] not in self.sample_set:
                    self.sample_set.append(row[SAMPLE_ID])

        # if there are multiple probes that agree, thats fine
        # need to remove conflicting and duplicate rows
        block_simple = self._simplify_block(block)
        if not block_simple:
            raise ConverterError(
                f"Oversimplified block {block[0]['Chr']}:{block[0]['Position']}"
            )
        block = block_simple
        # can assume each row in block has:
        #  same chm, pos
        #  unique sample names
        #  had consistent calls

        snp_names = tuple(sorted(frozenset([r[SNP_NAME] for r in block])))

        chm = block[0]["Chr"]
        # force chr prefix
        if not chm.startswith("chr"):
            chm = f"chr{chm}"
        if chm not in dbsnp_38_chrs_inv:
            raise ConverterError(
                f"Unexpected chromosome {block[0]['Chr']}:{block[0]['Position']}"
            )

        pos = int(block[0]["Position"])

        ref, alt = self.ref_alt_lookup(chm, pos)
        if not ref or not alt:
            raise ConverterError(
                f"Not a dbSNP SNV {block[0]['Chr']}:{block[0]['Position']}"
            )

        # using alts from dbSNP lists extra alts that aren't in samples or on the chip
        # using alts from called probes lists too few alts that are on the chip but not seen

        probed = set()
        for row in block:
            # get probe options & split
            # [A/C]
            probes = list(row["SNP"][1:-1].split("/"))
            # exclude indel probes
            if probes[0] not in "ATCG" or probes[1] not in "ATCG":
                raise ConverterError(
                    f"Not a SNV probe {block[0]['Chr']}:{block[0]['Position']}"
                )

            probed.update(probes)

        # strand correction if necessary ?
        # TODO determine this better!
        if ref not in probed:
            probed = set((strandswap[p] for p in probed))

        # hande ref/alt split
        if ref not in probed:
            raise ConverterError(
                f"Reference not probed {block[0]['Chr']}:{block[0]['Position']}"
            )
        probed.remove(ref)
        alt = tuple(sorted(probed))
        # alt may not be used, but is what the microarray could check for

        # convert calls
        calls = {}
        for row in block:

            allele1 = row[ALLELE1_FORWARD]
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

            allele2 = row[ALLELE2]
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

            assert row[SAMPLE_ID] not in calls
            calls[row[SAMPLE_ID]] = f"{allele1n}/{allele2n}"

        # always need to have all samples on all rows
        # even if that samples has been filtered out e.g. conflicting probes
        samples = []
        for sampleid in self.sample_set:
            samples.append({"GT": calls.get(sampleid, "./.")})

        vcfline = VCFLine(
            "", "", "", {}, chm, pos, snp_names, ref, alt, ".", ["PASS"], {}, samples
        )
        vcfline = self.clean_vcf_line(vcfline)

        return vcfline

    def clean_vcf_line(self, vcfline: VCFLine) -> VCFLine:
        # if no samples are called, discard
        if frozenset((s.get("GT", "./.") for s in vcfline.sample)) == frozenset(
            ("./.",)
        ):
            raise ConverterError(f"No calls made {vcfline.chrom}:{vcfline.pos}")

        return vcfline

    def convert(self, input, outfile):
        for vcfline in self._generate_vcf_lines(input):
            outfile.write(str(vcfline))
            outfile.write("\n")
            # outfile.flush()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dbsnp", help="path to DBSNP VCF")
    args = parser.parse_args()

    converter = Converter(args.dbsnp, args.dbsnp + ".tbi")

    # read from stdin
    # write to stdout
    converter.convert(sys.stdin, sys.stdout)
