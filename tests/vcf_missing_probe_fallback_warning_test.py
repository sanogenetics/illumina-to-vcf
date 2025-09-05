import io
import logging

import pytest

from illumina2vcf import IlluminaReader, VCFMaker
from illumina2vcf.bpm.bpmreader import CSVManifestReader, ManifestFilter
from illumina2vcf.bpm.referencegenome import ReferenceGenome


def _with_alias_snp_name(tab_report: io.StringIO, alias_suffix: str = "-alias") -> io.StringIO:
    """
    Copy a FinalReport and alter the first data row's SNP Name to an alias
    that won't match any manifest Name, to trigger the fallback path.
    """
    tab_report.seek(0)
    lines = tab_report.read().splitlines()
    out = []
    in_data = False
    header_idx = {}
    for i, line in enumerate(lines):
        out.append(line)
        if line == "[Data]":
            in_data = True
            continue
        if in_data and not header_idx:
            headers = line.split("\t")
            header_idx = {h: j for j, h in enumerate(headers)}
            continue
        if in_data and header_idx:
            parts = line.split("\t")
            if "SNP Name" in header_idx:
                parts[header_idx["SNP Name"]] = parts[header_idx["SNP Name"]] + alias_suffix
            out[-1] = "\t".join(parts)
            out.extend(lines[i + 1 :])
            break
    return io.StringIO("\n".join(out))


def test_missing_probe_falls_back_and_logs_warning(caplog):
    # Reference genome (trimmed test ref)
    genome_reader = ReferenceGenome(
        "tests/data/Homo_sapiens_assembly38.trim.fasta",
        "tests/data/Homo_sapiens_assembly38.trim.fasta.fai",
    )

    # Manifest (trimmed test manifest)
    bpm_records = ManifestFilter(frozenset(), skip_snps=False).filtered_records(
        CSVManifestReader(open("tests/data/GSA-24v3-0_A2.trim.csv"), genome_reader)
    )

    # Build a standard FinalReport then mutate the first SNP Name to an alias
    from .conftest import IlluminaBuilder

    report = IlluminaBuilder().build_file()
    alias_suffix = "-ALIAS123"
    report_alias = _with_alias_snp_name(report, alias_suffix=alias_suffix)

    # Parse and generate
    reader = IlluminaReader("\t")
    _ = reader.parse_header(report_alias)
    blocks = tuple(reader.generate_line_blocks(report_alias))

    vcfgenerator = VCFMaker(genome_reader, bpm_records)

    # Capture warnings from the module logger
    caplog.set_level(logging.WARNING, logger="illumina2vcf.vcf")

    # Generate some lines; ensure no crash and at least one non-comment VCF record is produced
    produced = 0
    for vline in vcfgenerator.generate_lines(blocks):
        if vline.comment_key or vline.comment_raw:
            continue
        produced += 1
        if produced > 20:
            break

    assert produced > 0

    # Confirm a warning was logged about falling back due to missing manifest probe name
    messages = "\n".join(m for _, _, m in caplog.record_tuples)
    assert "falling back to row-derived alleles" in messages
    assert "Probe name" in messages

