"""Conversion and release-safety tests, using synthetic participant identifiers only.

These tests use real VCFMaker logic, with a tiny fake genome supplying reference
alleles. Tests marked xfail are unmet safety requirements, not passing controls.
"""
from types import SimpleNamespace

import pytest

from illumina2vcf.illumina import IlluminaRow
from illumina2vcf.vcf import ConverterError, VCFMaker


class Genome:
    def __init__(self, bases=None):
        self.bases = bases or {("1", 100): "A", ("1", 200): "C", ("1", 300): "A"}
        self.reference_fasta = SimpleNamespace(
            faidx=SimpleNamespace(index={"1": SimpleNamespace(rlen=100000)}))

    def get_reference_bases(self, chrom, start, end):
        return self.bases.get((chrom, start + 1), "") if end == start + 1 else ""


def row(sample="P1", pos=100, snp="[A/G]", a1="A", a2="G", probe="rs1", chrom="1", strand="+"):
    return IlluminaRow(sample, probe, chrom, pos, snp, a1, a2, strand)


def records(maker, blocks):
    return [line for line in maker.generate_lines(blocks) if not (line.comment_key or line.comment_raw)]


@pytest.mark.parametrize("a1,a2,expected", [
    ("A", "A", "0/0"), ("A", "G", "0/1"), ("G", "G", "1/1"), ("-", "-", "./.")
])
def test_snp_genotype_is_preserved(a1, a2, expected):
    maker = VCFMaker(Genome())
    result = maker._line_block_to_vcf_line([row(a1=a1, a2=a2)], ["P1"])
    fields = str(result).split("\t")

    assert result.sample[0]["GT"] == expected
    assert result.ref == "A"
    assert tuple(result.alt) == ("G",)
    assert fields[6] == ("FAIL" if expected == "./." else "PASS")


def test_multiple_samples_do_not_cross_contaminate():
    maker = VCFMaker(Genome())
    result = maker._line_block_to_vcf_line([
        row(sample="P1", a1="G", a2="G"), row(sample="P2", a1="A", a2="A")], ["P1", "P2"])
    assert [s["GT"] for s in result.sample] == ["1/1", "0/0"]


def test_conflicting_identical_probe_calls_are_not_invented():
    maker = VCFMaker(Genome())
    result = maker._line_block_to_vcf_line([
        row(a1="A", a2="A"), row(a1="G", a2="G")], ["P1"])
    assert result.sample[0]["GT"] == "./."


def test_unprobed_reference_is_rejected():
    maker = VCFMaker(Genome())
    with pytest.raises(ConverterError, match="Reference .* not probed"):
        maker._line_block_to_vcf_line([row(snp="[C/G]", a1="C", a2="G")], ["P1"])


def test_invalid_chromosome_is_rejected():
    maker = VCFMaker(Genome())
    with pytest.raises(ConverterError):
        maker._line_block_to_vcf_line([row(chrom="unknown")], ["P1"])


def test_missing_indel_manifest_is_rejected():
    maker = VCFMaker(Genome({("1", 100): "A"}), bpm_records={})
    with pytest.raises(ConverterError, match="No BPM record for indel"):
        maker._line_block_to_vcf_line([row(snp="[I/D]", a1="I", a2="D")], ["P1"])


def test_valid_loci_continue_after_invalid_locus(caplog):
    maker = VCFMaker(Genome())
    result = records(maker, [[row(pos=100)], [row(pos=200, snp="[A/G]")],
                             [row(pos=300)]])
    assert [v.pos for v in result] == [100, 300]
    assert "Reference" in caplog.text



@pytest.mark.xfail(strict=True, reason="Current VCF header uses sample IDs only from the first locus")
def test_late_sample_must_be_in_vcf_column_header():
    maker = VCFMaker(Genome())
    output = list(maker.generate_lines([
        [row(sample="P1", pos=100)],
        [row(sample="P1", pos=300), row(sample="P2", pos=300, a1="G", a2="G")],
    ]))
    header = next(str(line) for line in output if line.comment_raw and "CHROM" in str(line))
    assert "P2" in header.split("\t")


@pytest.mark.parametrize("alleles,is_deletion,expected", [
    (("I", "I"), False, "1/1"), (("D", "D"), False, "0/0"),
    (("I", "D"), False, "0/1"), (("I", "I"), True, "0/0"),
    (("D", "D"), True, "1/1"), (("I", "D"), True, "0/1"),
])
def test_indel_genotype_encoding(alleles, is_deletion, expected):
    maker = VCFMaker(Genome())
    assert maker.convert_indel_genotype_to_vcf(alleles, is_deletion) == expected
