import csv
import logging
import re
from typing import Dict, Generator, Iterable, List, Tuple, Union

from fsspec.core import OpenFile
from puretabix.vcf import VCFLine
from pyfaidx import Fasta

from .illumina import ALLELE1, ALLELE2, SAMPLE_ID, SNP, SNP_NAME, STRAND
from .IlluminaBeadArrayFiles import RefStrand

STRANDSWAP = {"A": "T", "T": "A", "C": "G", "G": "C", "I": "I", "D": "D"}

logger = logging.getLogger(__name__)


class ConverterError(Exception):
    pass


class VCFMaker:

    reference: Union[str, OpenFile]
    reference_index: Union[str, OpenFile]
    reference_fasta: Fasta
    _buildsizes: Dict[str, str] = {}

    def __init__(self, genome_reader, indel_records) -> None:
        self._genome_reader = genome_reader
        self._indel_records = indel_records

    @staticmethod
    def is_valid_chromosome(chrom: str):
        return re.match(r"^chr[1-9XYM][0-9]?$", chrom)

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
        for chrom, rec in self._genome_reader.reference_fasta.faidx.index.items():
            if self.is_valid_chromosome(chrom):
                yield VCFLine.as_comment_key_dict(
                    "contig",
                    {"ID": chrom, "length": rec.rlen, "assembly": buildname},
                )

    def generate_lines(self, blocks) -> Generator[VCFLine, None, None]:
        column_header = False
        samples: List[str] = []
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
        # if block is silly big skip it (first block will be included regardless of size)
        if sample_set and len(block) > 10 * len(sample_set):
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
        if chm not in self._genome_reader.reference_fasta.faidx.index:
            raise ConverterError(f"Unexpected chromosome {chm}:{block[0]['Position']}")
        if not self.is_valid_chromosome(chm):
            raise ConverterError(f"Unexpected chromosome {chm}:{block[0]['Position']}")
        pos = int(block[0]["Position"])

        ref = self._genome_reader.get_reference_bases(chm, pos - 1, pos)
        converted_calls = {}
        # handel indels
        if "I" in probed or "D" in probed:
            if probed.intersection({"A", "C", "G", "T"}):
                raise ConverterError(
                    f"{';'.join(snp_names)}: contains indels and SNPs ({','.join(probed)}) {block[0]['Chr']}:{block[0]['Position']}"
                )
            try:
                locus_records = self._indel_records[chm][pos]
                (ref, alt) = self.get_alleles_for_indel(locus_records[0])
                for alt_record in locus_records[1:]:
                    if (ref, alt) != self.get_alleles_for_indel(alt_record):
                        raise ConverterError(
                            f"{';'.join(snp_names)}: Mismatched alleles ({','.join(probed)}) {block[0]['Chr']}:{block[0]['Position']}"
                        )
                # Illumina position is the position of the insertion (ie; one base after the ref allele)
                # We want the position of the start of the ref allele, so we need to reduce pos by 1
                pos -= 1

            except KeyError:
                raise ConverterError(
                    f"{';'.join(snp_names)}: No BPM record ({','.join(probed)}) {block[0]['Chr']}:{block[0]['Position']}"
                )
            for sampleid in calls:
                converted_calls[sampleid] = self.convert_indel_genotype_to_vcf(
                    calls[sampleid], locus_records[0].is_deletion
                )
        else:
            # handle ref/alt split
            if ref not in probed:
                raise ConverterError(
                    f"{';'.join(snp_names)}: Reference ({ref}) not probed ({','.join(probed)}) {block[0]['Chr']}:{block[0]['Position']}"
                )
            probed.remove(ref)

            # alt may not be used, but is what the microarray could check for
            alt = tuple(sorted(probed))

            # convert calls
            for sampleid in calls:
                allele1 = calls[sampleid][0]
                if allele1 not in "ATCG":
                    allele1n = "."
                elif allele1 == ref:
                    allele1n = "0"
                elif allele1 not in alt:
                    raise ConverterError(f"Unexpected forward {block[0]['Chr']}:{block[0]['Position']}")
                else:
                    allele1n = str(1 + alt.index(allele1))

                allele2 = calls[sampleid][1]
                if allele2 not in "ATCG":
                    allele2n = "."
                elif allele2 == ref:
                    allele2n = "0"
                elif allele2 not in alt:
                    raise ConverterError(f"Unexpected forward {block[0]['Chr']}:{block[0]['Position']}")
                else:
                    allele2n = str(1 + alt.index(allele2))

                assert sampleid not in converted_calls
                converted_calls[sampleid] = self.format_vcf_genotype(allele1n, allele2n)

        # always need to have all samples on all rows
        # even if that samples has been filtered out e.g. conflicting probes
        samples = []
        for sampleid in sample_set:
            samples.append({"GT": converted_calls.get(sampleid, "./.")})

        vcfline = VCFLine("", "", "", {}, chm, pos, snp_names, ref, alt, ".", ["PASS"], {}, samples)

        return vcfline

    def get_alleles_for_indel(self, bpm_record):
        """
            Determine REF and ALT alleles for indel manifest record
            Args:
               bpm_record: BPM record for indel
        Returns:
            tupple of REF and ALT alleles
        """

        (_, indel_sequence, _) = bpm_record.get_indel_source_sequences(RefStrand.Plus)
        start_index = bpm_record.pos - 1
        chrom = bpm_record.chromosome
        if chrom == "XX" or chrom == "XY":
            chrom = "X"

        if bpm_record.is_deletion:
            reference_base = self._genome_reader.get_reference_bases(chrom, start_index - 1, start_index)
            reference_allele = reference_base + indel_sequence
            alternate_allele = reference_base
        else:
            reference_base = self._genome_reader.get_reference_bases(chrom, start_index, start_index + 1)
            reference_allele = reference_base
            alternate_allele = reference_base + indel_sequence
        return (reference_allele, alternate_allele)

    def format_vcf_genotype(self, vcf_allele1_char, vcf_allele2_char):
        """
        Create a VCF representation of the genotype based on alleles. Format appropriately if haploid.
        Args:
            vcf_allele1_char (string): 0,1,2 etc.
            vcf_allele2_char (string): 0,1,2, etc.
        Returns
            string: String representation of genotype (e.g., "0/1")
        """

        vcf_genotype = ""
        if vcf_allele2_char < vcf_allele1_char:
            vcf_genotype = str(vcf_allele2_char) + "/" + str(vcf_allele1_char)
        else:
            vcf_genotype = str(vcf_allele1_char) + "/" + str(vcf_allele2_char)
        return vcf_genotype

    def convert_indel_genotype_to_vcf(self, nucleotide_genotypes, is_deletion):
        """
        For indel, convert indel SNP genotype (e.g., I/D) into VCF genotype (e.g, 0/1)
        Args:
            nucleotide_genotypes (string,string): SNP genotype from manifest (e.g., ('D', 'I'))
            is_deletion (bool): Whether the BPM record that produced the nucleotide genotypes is a reference deletion
        Returns:
            string: VCF genotype (e.g, "0/1")
        """
        if not nucleotide_genotypes:
            return self.format_vcf_genotype(".", ".")

        if is_deletion:
            vcf_allele1_char = "0" if nucleotide_genotypes[0] == "I" else "1"
            vcf_allele2_char = "0" if nucleotide_genotypes[1] == "I" else "1"
        else:
            vcf_allele1_char = "1" if nucleotide_genotypes[0] == "I" else "0"
            vcf_allele2_char = "1" if nucleotide_genotypes[1] == "I" else "0"

        vcf_genotype = self.format_vcf_genotype(vcf_allele1_char, vcf_allele2_char)

        return vcf_genotype

    def _simplify_block(self, block: List[Dict[str, str]]) -> Tuple[Dict[str, Tuple[str, str]], set]:
        combined_probes = {}
        calls = {}
        conflicts = []
        for row in block:
            sampleid = row[SAMPLE_ID]
            strand = row[STRAND]

            probes = row[SNP]
            match = re.match(r"\[([ATCGID]+)\/([ATCGID]+)\]", probes)
            if not match:
                raise ConverterError(f"Unexpcted probes {row['Chr']}:{row['Position']} {probes}")

            new_probes = match.group(1, 2)
            # exclude indel probes
            # if new_probes[0] not in "ATCG" or new_probes[1] not in "ATCG":
            #    raise ConverterError(f"Not a SNV probe {row['Chr']}:{row['Position']} {probes}")

            # swap strand if needed
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
                        calls[sampleid], new_calls, tuple(combined_probes[sampleid]), new_probes
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
        previous_probes: Tuple[str, ...],
        new_probes: Tuple[str, ...],
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
