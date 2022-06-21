import csv
import logging
from typing import Dict, Generator, Iterable, List, Tuple

from puretabix.vcf import VCFLine
from pyfaidx import Fasta

from .illumina import ALLELE1, ALLELE2, SAMPLE_ID, SNP_NAME, STRAND

STRANDSWAP = {"A": "T", "T": "A", "C": "G", "G": "C"}

logger = logging.getLogger(__name__)


class ConverterError(Exception):
    pass


class VCFMaker:

    reference_filename: str
    reference_index_filename: str
    reference: Fasta
    _buildsizes: Dict[str, str] = {}

    def __init__(self, reference_filename: str, reference_index_filename: str) -> None:
        self.reference_filename = reference_filename
        self.reference_index_filename = reference_index_filename
        self.reference = Fasta(reference_filename)

    @property
    def buildsizes(self) -> Dict[str, str]:
        if not self._buildsizes:
            with open(self.reference_index_filename, "rt") as index_fh:
                index_reader = csv.reader(index_fh, delimiter="\t")
                for line in index_reader:
                    self._buildsizes[line[0]] = line[1]
        return self._buildsizes

    def ref_lookup(self, chm: str, pos: int) -> str:
        # VCF is 1-based
        # FASTA is 0-based
        ref_base = str(self.reference[chm][pos - 1])
        return ref_base

    def generate_header(self, date: str, source: str, buildname: str) -> Generator[VCFLine, None, None]:
        # write header
        yield VCFLine.as_comment_key_string("fileformat", "VCFv4.3")
        yield VCFLine.as_comment_key_string("filedate", date)
        yield VCFLine.as_comment_key_string("source", f'"{source}, Sano Genetics"')

        # ##FILTER=<ID=PASS,Description="All filters passed">
        yield VCFLine.as_comment_key_dict("FILTER", {"ID": "PASS", "Description": '"All filters passed"'})

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

        # ##contig=<ID=1,length=249250621,assembly=GRCh37>
        for chrom, length in self.buildsizes.items():
            yield VCFLine.as_comment_key_dict(
                "contig",
                {"ID": chrom, "length": length, "assembly": buildname},
            )

    def generate_lines(self, blocks) -> Generator[VCFLine, None, None]:
        column_header = False
        samples = []
        for block in blocks:
            # if we've not got a list of samples yet, get them from this block unfiltered
            if not samples:
                for row in block:
                    if row[SAMPLE_ID] not in samples:
                        samples.append(row[SAMPLE_ID])
            # if we've not output column names, do so
            if not column_header:
                yield self._get_column_header(samples)
                column_header = True

            try:
                yield self._line_block_to_vcf_line(block, samples)
            except ConverterError as e:
                logger.error(e)
                continue

    def _get_column_header(self, samples: Iterable[str]) -> VCFLine:
        return VCFLine.as_comment_raw(
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
                    "\t".join(samples),  # customizable name of sample
                )
            )
        )

    def _line_block_to_vcf_line(self, block: List[Dict[str, str]], sample_set) -> VCFLine:
        # if block is silly big skip it (this will fail with >100 samples)
        if len(block) > 100 or (sample_set and len(block) > 10 * len(sample_set)):
            raise ConverterError(f"Oversized block {block[0]['Chr']}:{block[0]['Position']}: {len(block)} rows")

        # if we've not got a list of samples yet, get them from this block unfiltered
        if not sample_set:
            for row in block:
                if row[SAMPLE_ID] not in sample_set:
                    sample_set.append(row[SAMPLE_ID])

        # if there are multiple probes that agree, thats fine
        # need to remove conflicting and duplicate rows
        [calls, probed] = self._simplify_block(block)
        # if not calls:
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
        if chm not in self.buildsizes.keys():
            raise ConverterError(f"Unexpected chromosome {block[0]['Chr']}:{block[0]['Position']}")
        pos = int(block[0]["Position"])

        ref = self.ref_lookup(chm, pos)

        # hande ref/alt split
        if ref not in probed:
            raise ConverterError(
                f"{';'.join(snp_names)}: Reference ({ref}) not probed ({','.join(probed)}) {block[0]['Chr']}:{block[0]['Position']}"
            )
        probed.remove(ref)

        # alt may not be used, but is what the microarray could check for
        alt = tuple(sorted(probed))

        # convert calls
        converted_calls = {}
        for sampleid in calls:
            allele1 = calls[sampleid][0]
            if allele1 not in "ATCG":
                allele1n = "."
            elif allele1 == ref:
                allele1n = 0
            elif allele1 not in alt:
                raise ConverterError(f"Unexpected forward {block[0]['Chr']}:{block[0]['Position']}")
            else:
                allele1n = 1 + alt.index(allele1)

            allele2 = calls[sampleid][1]
            if allele2 not in "ATCG":
                allele2n = "."
            elif allele2 == ref:
                allele2n = 0
            elif allele2 not in alt:
                raise ConverterError(f"Unexpected forward {block[0]['Chr']}:{block[0]['Position']}")
            else:
                allele2n = 1 + alt.index(allele2)

            assert sampleid not in converted_calls
            converted_calls[sampleid] = f"{allele1n}/{allele2n}"

        # always need to have all samples on all rows
        # even if that samples has been filtered out e.g. conflicting probes
        samples = []
        for sampleid in sample_set:
            samples.append({"GT": converted_calls.get(sampleid, "./.")})

        vcfline = VCFLine("", "", "", {}, chm, pos, snp_names, ref, alt, ".", ["PASS"], {}, samples)
        # vcfline = self.clean_vcf_line(vcfline)

        return vcfline

    def _simplify_block(self, block: List[Dict[str, str]]) -> Tuple[Dict[str, Tuple[str, str]], set]:
        combined_probes = {}
        calls = {}
        conflicts = []
        for row in block:
            sampleid = row[SAMPLE_ID]
            strand = row[STRAND]

            new_probes = tuple(row["SNP"][1:-1].split("/"))
            # exclude indel probes
            if new_probes[0] not in "ATCG" or new_probes[1] not in "ATCG":
                raise ConverterError(f"Not a SNV probe {row['Chr']}:{row['Position']}")
            if strand == "-":
                new_probes = tuple((STRANDSWAP[probe] for probe in new_probes))
            new_calls = (row[ALLELE1], row[ALLELE2])

            if sampleid not in calls:
                # This sample hasnt been called before
                calls[sampleid] = new_calls
                combined_probes[sampleid] = set(new_probes)
            else:
                # This sample has been called before
                # need to reconcile multiple calls
                try:
                    calls[sampleid] = self._combine_calls(
                        calls[sampleid], new_calls, combined_probes[sampleid], new_probes
                    )
                except ConverterError as e:
                    conflicts.append(sampleid)
                combined_probes[sampleid].update(new_probes)

        probed = set().union(*[combined_probes[sampleid] for sampleid in combined_probes.keys()])

        # for all samples with conflicting genotype calls, set the call to no call
        for sampleid in conflicts:
            calls[sampleid] = ("-", "-")

        return calls, probed

    def _combine_calls(
        self,
        previous_calls: Tuple[str, str],
        new_calls: Tuple[str, str],
        previous_probes: Tuple[str, str],
        new_probes: Tuple[str, str],
    ):
        if previous_calls == new_calls or previous_calls == new_calls[::-1] or new_calls == ("-", "-"):
            return previous_calls
        elif new_calls[0] == new_calls[1] and new_calls[0] in previous_calls:
            if not new_probes[0] in previous_calls or not new_probes[1] in previous_calls:
                return previous_calls  # conflicting call is a homozygote where only one of the probes matches the true genotype
            else:
                raise ConverterError("conflict!")  # both probes match but different call
        elif previous_calls == ("-", "-"):
            return new_calls
        elif previous_calls[0] == previous_calls[1] and previous_calls[0] in new_calls:
            if new_probes[0] not in previous_probes or new_probes[1] not in previous_probes:
                return new_calls  # conflicting call is a het where only one of the probes has been tested previously (and previous call is homozygous)
            else:
                raise ConverterError("conflict!")  # both probes match but different call
        else:
            raise ConverterError("conflict!")  # conflicting hets or conflicting homozygotes