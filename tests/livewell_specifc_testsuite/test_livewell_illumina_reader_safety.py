"""Livewell input integrity checks for illumina2vcf. Run from repository root."""
from dataclasses import replace
from io import StringIO

import pytest

from illumina2vcf.illumina import DataError, IlluminaReader

COLUMNS = ["Sample ID", "SNP Name", "Chr", "Position", "SNP", "Allele1 - Plus", "Allele2 - Plus", "Plus/Minus Strand"]


def report(*records):
    return StringIO("\t".join(COLUMNS) + "\n" + "\n".join("\t".join(map(str, r)) for r in records) + "\n")


def row(sample="P1", probe="rs1", chrom="1", pos=100, snp="[A/G]", a1="A", a2="G", strand="+"):
    return (sample, probe, chrom, pos, snp, a1, a2, strand)


def test_input_reordering_does_not_change_locus_groups():
    source = [row("P2", "rs2", pos=200), row("P1", "rs1"), row("P2", "rs1"), row("P1", "rs2", pos=200)]
    reader = IlluminaReader("\t")
    def canonical(rows):
        return [[(r.sample_id, r.snp_name, r.chrom, r.pos, r.allele1, r.allele2) for r in block]
                for block in reader.generate_line_blocks(report(*rows))]
    assert canonical(source) == canonical(list(reversed(source)))


def test_blocklist_removes_only_target_probe(tmp_path):
    blocked = tmp_path / "blocklist.txt"
    blocked.write_text("bad_probe\n", encoding="utf-8")
    blocks = list(IlluminaReader("\t", str(blocked)).generate_line_blocks(
        report(row(probe="bad_probe"), row(probe="good_probe"))))
    assert [[r.snp_name for r in b] for b in blocks] == [["good_probe"]]


def test_illumina_duplicates_and_zero_position_are_excluded():
    blocks = list(IlluminaReader("\t").generate_line_blocks(report(
        row(probe="rs1ilmndup"), row(probe="zero", pos=0), row(probe="valid"))))
    assert [[r.snp_name for r in b] for b in blocks] == [["valid"]]


@pytest.mark.parametrize("missing", ["Chr", "Position", "SNP Name", "SNP"])
def test_missing_required_column_is_rejected(missing):
    columns = [c for c in COLUMNS if c != missing]
    source = StringIO("\t".join(columns) + "\n" + "\t".join(
        str(value) for col, value in zip(COLUMNS, row()) if col != missing) + "\n")
    with pytest.raises(DataError):
        list(IlluminaReader("\t").generate_line_blocks(source))


def test_nonnumeric_position_is_rejected():
    with pytest.raises(ValueError):
        list(IlluminaReader("\t").generate_line_blocks(report(row(pos="unknown"))))


@pytest.mark.xfail(strict=True, reason="Safety requirement: duplicate sample/probe should be explicitly rejected; current reader retains duplicates")
def test_duplicate_sample_probe_must_not_be_silently_accepted():
    source = report(row(), row())
    with pytest.raises(DataError):
        list(IlluminaReader("\t").generate_line_blocks(source))


@pytest.mark.xfail(strict=True, reason="Safety requirement: missing Sample ID must fail; current reader does not validate blank values")
def test_blank_sample_id_is_rejected():
    with pytest.raises(DataError):
        list(IlluminaReader("\t").generate_line_blocks(report(row(sample=""))))
