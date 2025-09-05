import re
import pytest

from illumina2vcf import IlluminaReader, VCFMaker
from illumina2vcf.bpm.referencegenome import ReferenceGenome
from illumina2vcf.vcf import ConverterError


def _is_indel_block(block) -> bool:
    """Return True if the block's alleles are indel-only (I/D)."""
    alleles = set()
    for row in block:
        m = re.match(r"\[([ATCGID]+)/(\w+)\]", row.snp)
        if not m:
            continue
        a1, a2 = m.group(1), m.group(2)
        alleles.update([a1, a2])
    # indel-only if all alleles are subset of {I, D}
    return bool(alleles) and alleles.issubset({"I", "D"})


def test_indel_block_without_bpm_indel_record_raises():
    # Build reference
    genome_reader = ReferenceGenome(
        "tests/data/Homo_sapiens_assembly38.trim.fasta",
        "tests/data/Homo_sapiens_assembly38.trim.fasta.fai",
    )

    # Create blocks from a synthetic report
    from .conftest import IlluminaBuilder

    reader = IlluminaReader("\t")
    report = IlluminaBuilder().build_file()
    _ = reader.parse_header(report)
    blocks = tuple(reader.generate_line_blocks(report))

    # Pick an indel-only block from the synthetic data
    indel_blocks = [b for b in blocks if _is_indel_block(b)]
    assert indel_blocks, "Expected at least one indel-only block in test data"
    block = indel_blocks[0]

    # Use an empty bpm_records dict so no BPM record is available for the indel locus
    vcfgenerator = VCFMaker(genome_reader, bpm_records={})

    # Expect a clear failure when trying to create a VCF line for an indel without BPM context
    with pytest.raises(ConverterError) as ei:
        # Call the internal to hit the specific error deterministically
        vcfgenerator._line_block_to_vcf_line(block, sample_set=[])

    assert "No BPM record for indel" in str(ei.value)

