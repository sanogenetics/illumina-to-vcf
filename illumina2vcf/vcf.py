import contextlib
import logging
import re
from typing import Dict, Generator, Iterable, List, Literal, Optional, Tuple

from puretabix.vcf import VCFLine

from illumina2vcf.bpm.bpmrecord import BPMRecord
from illumina2vcf.bpm.illuminabeadarrayfiles import RefStrand
from illumina2vcf.bpm.referencegenome import ReferenceGenome
from illumina2vcf.illumina import IlluminaRow

STRANDSWAP = {"A": "T", "T": "A", "C": "G", "G": "C", "I": "I", "D": "D"}

logger = logging.getLogger(__name__)


class ConverterError(Exception):
    pass


class QcStats:
    def __init__(self):
        self.vcf_lines_autosomal = 0
        self.called_lines_autosomal = 0
        self.heterozygous_lines_autosomal = 0
        self.called_lines_x = 0
        self.heterozygous_lines_x = 0
        self.vcf_lines_y = 0
        self.called_lines_y = 0

    def vcf_comment(self) -> Optional[VCFLine]:
        qc_stats = {}
        with contextlib.suppress(ZeroDivisionError):
            qc_stats["callrate"] = f"{self.call_rate('autosomal'):.2f}"
        with contextlib.suppress(ZeroDivisionError):
            qc_stats["y_notnull"] = f"{self.call_rate('Y'):.2f}"
        with contextlib.suppress(ZeroDivisionError):
            qc_stats["het"] = f"{self.heterozygosity('autosomal'):.2f}"
        with contextlib.suppress(ZeroDivisionError):
            qc_stats["x_het"] = f"{self.heterozygosity('X'):.2f}"
        return VCFLine.as_comment_key_dict("qc_stats", qc_stats) if qc_stats else None

    def call_rate(self, chrom: Literal["autosomal", "Y"]) -> float:
        if chrom == "autosomal":
            return self.called_lines_autosomal / self.vcf_lines_autosomal
        elif chrom == "Y":
            return self.called_lines_y / self.vcf_lines_y
        else:
            msg = f"cannot calculate call rate for chrom {chrom}"
            raise NotImplementedError(msg)

    def heterozygosity(self, chrom: Literal["autosomal", "X"]) -> float:
        if chrom == "autosomal":
            return self.heterozygous_lines_autosomal / self.called_lines_autosomal
        elif chrom == "X":
            return self.heterozygous_lines_x / self.called_lines_x
        else:
            msg = f"cannot calculate heterozygosity for chrom {chrom}"
            raise NotImplementedError(msg)

    def _update_qc_stats(self, vcfline: VCFLine) -> None:
        gt = vcfline.sample[0]["GT"]
        chm = vcfline.chrom
        if chm == "chrM":
            pass
        elif chm == "chrX":
            if gt != "./.":
                self.called_lines_x += 1
                if len(set(gt.split("/"))) == 2:  # noqa: PLR2004
                    self.heterozygous_lines_x += 1
        elif chm == "chrY":
            self.vcf_lines_y += 1
            if gt != "./.":
                self.called_lines_y += 1
        else:
            self.vcf_lines_autosomal += 1
            if gt != "./.":
                self.called_lines_autosomal += 1
                if len(set(gt.split("/"))) == 2:  # noqa: PLR2004
                    self.heterozygous_lines_autosomal += 1


class VCFMaker:
    _genome_reader: ReferenceGenome
    _indel_records: Dict

    def __init__(self, genome_reader: ReferenceGenome, indel_records: Optional[Dict] = None):
        self._indel_records = {} if indel_records is None else indel_records
        self._genome_reader = genome_reader
        self.qc_stats = QcStats()

    @staticmethod
    def is_valid_chromosome(chrom: str) -> bool:
        return bool(re.match(r"^chr[1-9XYM][0-9]?$", chrom))

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
        # ##qc_stats=<callrate=0.99,het=0.33,x_het=0.21,y_notnull=0.24>
        if qc_stats := self.qc_stats.vcf_comment():
            yield qc_stats

    def generate_lines(self, blocks: Iterable[List[IlluminaRow]]) -> Generator[VCFLine, None, None]:
        column_header = False
        samples: List[str] = []
        for block in blocks:
            # if we've not got a list of samples yet, get them from this block unfiltered
            if not samples:
                for row in block:
                    if row.sample_id not in samples:
                        samples.append(row.sample_id)
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

    def _line_block_to_vcf_line(self, block: List[IlluminaRow], sample_set) -> VCFLine:
        # if block is silly big skip it (first block will be included regardless of size)
        if sample_set and len(block) > 10 * len(sample_set):
            msg = f"Oversized block {block[0].chrom}:{block[0].pos}: {len(block)} rows"
            raise ConverterError(msg)

        # if we've not got a list of samples yet, get them from this block unfiltered
        if not sample_set:
            for row in block:
                if row.sample_id not in sample_set:
                    sample_set.append(row.sample_id)

        # if there are multiple probes that agree, thats fine
        # need to remove conflicting and duplicate rows
        [calls, probed] = self._simplify_block(block)
        # if not calls:
        #    raise ConverterError(
        #        f"Oversimplified block {block[0].chrom}:{block[0].pos}"
        #    )
        # can assume each row in block has:
        #  same chm, pos
        #  unique sample names
        #  had consistent calls

        snp_names = tuple(sorted(frozenset(r.snp_name for r in block)))

        chm = block[0].chrom
        # convert pseudoautosomal (XY) to X
        if chm in ("XY", "chrXY"):
            chm = "chrX"
        # convert MT to M
        if chm in ("MT", "chrMT"):
            chm = "chrM"
        # force chr prefix
        if not chm.startswith("chr"):
            chm = f"chr{chm}"
        if chm not in self._genome_reader.reference_fasta.faidx.index:
            msg = f"Unexpected chromosome {chm}:{block[0].pos}"
            raise ConverterError(msg)
        if not self.is_valid_chromosome(chm):
            msg = f"Invalid chromosome {chm}:{block[0].pos}"
            raise ConverterError(msg)
        pos = int(block[0].pos)

        ref = self._genome_reader.get_reference_bases(chm, pos - 1, pos)
        if not ref:
            msg = f"Unable to get reference {chm}:{pos-1}-{pos}"
            raise ConverterError(msg)

        converted_calls = {}
        # handel indels
        if "I" in probed or "D" in probed:
            if probed.intersection({"A", "C", "G", "T"}):
                msg = f"{';'.join(snp_names)}: contains indels and SNPs ({','.join(probed)}) {block[0].chrom}:{block[0].pos}"
                raise ConverterError(msg)

            if (chm, pos) in self._indel_records:
                # get the records in the manifest for this locaion
                locus_records = self._indel_records[(chm, pos)]

                (ref, alt, pos) = self.get_alleles_for_indel(locus_records[0])

                # check the other lotus records have the same reference + alternative alleles
                for alt_record in locus_records[1:]:
                    if (ref, alt, pos) != self.get_alleles_for_indel(alt_record):
                        msg = f"{';'.join(snp_names)}: Mismatched alleles ({','.join(probed)}) {block[0].chrom}:{block[0].pos}"
                        raise ConverterError(msg)
                # Illumina position is the position of the insertion (ie; one base after the ref allele)
                # We want the position of the start of the ref allele, so we need to reduce pos by 1

                # alt column in VCF is a list
                alt = (alt,)

            else:
                msg = f"{';'.join(snp_names)}: No BPM record for indel ({','.join(probed)}) {chm}:{pos}"
                raise ConverterError(msg)

            for sampleid in calls:
                converted_calls[sampleid] = self.convert_indel_genotype_to_vcf(
                    calls[sampleid], locus_records[0].is_deletion
                )
        else:
            # SNPs
            # handle ref/alt split
            if ref not in probed:
                msg = f"{';'.join(snp_names)}: Reference ({ref}) not probed ({','.join(probed)}) {block[0].chrom}:{block[0].pos}"
                raise ConverterError(msg)
            probed.remove(ref)

            # alt may not be used, but is what the microarray could check for
            # alt column in VCF is a list
            alt = tuple(sorted(probed))

            # convert calls
            for sampleid in calls:
                allele1 = calls[sampleid][0]
                if allele1 not in "ATCG":
                    allele1n = "."
                elif allele1 == ref:
                    allele1n = "0"
                elif allele1 not in alt:
                    msg = f"Unexpected forward {block[0].chrom}:{block[0].pos} {allele1} vs {ref}/{alt}"
                    raise ConverterError(msg)
                else:
                    allele1n = str(1 + alt.index(allele1))

                allele2 = calls[sampleid][1]
                if allele2 not in "ATCG":
                    allele2n = "."
                elif allele2 == ref:
                    allele2n = "0"
                elif allele2 not in alt:
                    msg = f"Unexpected forward {block[0].chrom}:{block[0].pos} {allele2} vs {ref}/{alt}"
                    raise ConverterError(msg)
                else:
                    allele2n = str(1 + alt.index(allele2))

                assert sampleid not in converted_calls  # noqa: S101
                converted_calls[sampleid] = self.format_vcf_genotype(allele1n, allele2n)

        # always need to have all samples on all rows
        # even if that samples has been filtered out e.g. conflicting probes
        samples = [{"GT": converted_calls.get(sampleid, "./.")} for sampleid in sample_set]

        vcfline = VCFLine("", "", "", {}, chm, pos, snp_names, ref, alt, ".", ["PASS"], {}, samples)

        # update qc totals if single sample
        if len(samples) == 1:
            self.qc_stats._update_qc_stats(vcfline)

        return vcfline

    def get_alleles_for_indel(self, bpm_record: BPMRecord):
        """
            Determine REF and ALT alleles for indel manifest record
            Args:
               bpm_record: BPM record for indel
        Returns:
            tuple of REF and ALT alleles
        """

        (_, indel_sequence, _) = bpm_record.get_indel_source_sequences(RefStrand.Plus)
        start_index = bpm_record.pos - 1
        chrom = bpm_record.chromosome
        if chrom in ("XX", "XY"):
            chrom = "X"

        if bpm_record.is_deletion:
            reference_base = self._genome_reader.get_reference_bases(chrom, start_index - 1, start_index)
            reference_allele = reference_base + indel_sequence
            alternate_allele = reference_base
            pos = bpm_record.pos - 1
        else:
            reference_base = self._genome_reader.get_reference_bases(chrom, start_index, start_index + 1)
            reference_allele = reference_base
            alternate_allele = reference_base + indel_sequence
            pos = bpm_record.pos
        return (reference_allele, alternate_allele, pos)

    def format_vcf_genotype(self, vcf_allele1_char, vcf_allele2_char):
        """
        Create a VCF representation of the genotype based on alleles. Format appropriately if haploid.
        Args:
            vcf_allele1_char (string): 0,1,2 etc.
            vcf_allele2_char (string): 0,1,2, etc.
        Returns
            string: String representation of genotype (e.g., "0/1")
        """

        return (
            f"{vcf_allele2_char!s}/{vcf_allele1_char!s}"
            if vcf_allele2_char < vcf_allele1_char
            else f"{vcf_allele1_char!s}/{vcf_allele2_char!s}"
        )

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

        return self.format_vcf_genotype(vcf_allele1_char, vcf_allele2_char)

    def _simplify_block(self, block: List[IlluminaRow]) -> Tuple[Dict[str, Tuple[str, str]], set]:
        combined_probes = {}
        calls = {}
        conflicts = []
        for row in block:
            sampleid = row.sample_id
            strand = row.strand

            probes = row.snp
            match = re.match(r"\[([ATCGID]+)\/([ATCGID]+)\]", probes)
            if not match:
                msg = f"Unexpcted probes {row.chrom}:{row.pos} {probes}"
                raise ConverterError(msg)

            new_probes = match.group(1, 2)
            # exclude indel probes
            # if new_probes[0] not in "ATCG" or new_probes[1] not in "ATCG":
            #    raise ConverterError(f"Not a SNV probe {row.chrom}:{row.pos} {probes}")

            # swap strand if needed
            if strand == "-":
                new_probes = tuple(STRANDSWAP[probe] for probe in new_probes)

            new_calls = (row.allele1, row.allele2)

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
                except ConverterError:
                    conflicts.append(sampleid)
                combined_probes[sampleid].update(new_probes)

        probed = set().union(*[combined_probes[sampleid] for sampleid in combined_probes])

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
        if previous_calls in (new_calls, new_calls[::-1]):
            return previous_calls
        elif new_calls[0] == new_calls[1] and new_calls[0] in previous_calls:
            if new_probes[0] not in previous_calls or new_probes[1] not in previous_calls:
                return (
                    previous_calls
                )  # conflicting call is a homozygote where only one of the probes matches the true genotype
            else:
                raise ConverterError("conflict!")  # both probes match but different call  # noqa: EM101
        elif previous_calls == ("-", "-"):
            return new_calls
        elif previous_calls[0] == previous_calls[1] and previous_calls[0] in new_calls:
            if new_probes[0] not in previous_probes or new_probes[1] not in previous_probes:
                return new_calls  # conflicting call is a het where only one of the probes has been tested previously (and previous call is homozygous)
            else:
                raise ConverterError("conflict!")  # both probes match but different call  # noqa: EM101
        else:
            raise ConverterError("conflict!")  # conflicting hets or conflicting homozygotes  # noqa: EM101
