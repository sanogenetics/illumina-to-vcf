from typing import Dict, List, Tuple

from pytest import fixture

from illumina2vcf import IlluminaReader, VCFMaker
from illumina2vcf.bpm.BPMReader import CSVManifestReader, ManifestFilter
from illumina2vcf.bpm.BPMRecord import BPMRecord
from illumina2vcf.bpm.ReferenceGenome import ReferenceGenome

from .conftest import IlluminaBuilder


@fixture
def blocks() -> tuple[List[Dict[str, str]], ...]:
    reader = IlluminaReader("\t")
    illumina = IlluminaBuilder().build_file()
    _, _ = reader.parse_header(illumina)
    blocks = tuple(reader.generate_line_blocks(illumina))
    return blocks


@fixture
def genome_reader() -> ReferenceGenome:
    return ReferenceGenome(
        "tests/data/Homo_sapiens_assembly38.trim.fasta", "tests/data/Homo_sapiens_assembly38.trim.fasta.fai"
    )


@fixture
def indel_records(genome_reader) -> Dict[Tuple[str, int], BPMRecord]:
    return ManifestFilter(
        CSVManifestReader("tests/data/GSA-24v3-0_A2.trim.csv", genome_reader), frozenset(), skip_snps=True
    ).filtered_records()


class TestVCF:
    def test_header(self, genome_reader):
        """
        GIVEN a VCF generator
        """
        vcfgenerator = VCFMaker(genome_reader, {})
        """
        WHEN generating a header
        """
        lines = tuple((str(i) for i in vcfgenerator.generate_header("date", "source", "build")))
        """
        THEN it should produce valid header output
        """
        for line in lines:
            assert line
            assert line[0] == "#"

    def test_data(self, blocks, genome_reader):
        """
        GIVEN an interable of blocks
        """
        vcfgenerator = VCFMaker(genome_reader, {})
        """
        GIVEN a VCF generator
        """
        for _ in vcfgenerator.generate_header("date", "source", "build"):
            pass
        lines = vcfgenerator.generate_lines(blocks)
        """
        THEN it should produce valid output
        """
        maxref = 0
        maxalt = 0
        for line in lines:
            if line.comment_key or line.comment_raw:
                continue
            assert line.ref
            assert line.alt
            assert max((len(i) for i in line.alt))
            assert len(line.sample) == 1
            maxref = max(maxref, len(line.ref))
            maxalt = max(maxalt, max((len(i) for i in line.alt)))
        assert maxref == 1
        assert maxalt == 1

    def disabled_test_data_indel(self, blocks, genome_reader, indel_records):
        """
        GIVEN an interable of blocks and indel data
        """
        vcfgenerator = VCFMaker(genome_reader, indel_records)
        """
        GIVEN a VCF generator
        """
        for _ in vcfgenerator.generate_header("date", "source", "build"):
            pass
        lines = vcfgenerator.generate_lines(blocks)
        """
        THEN it should produce valid output
        """
        maxref = 0
        maxalt = 0
        for line in lines:
            if line.comment_key or line.comment_raw:
                continue
            assert line.ref
            assert line.alt
            assert max((len(i) for i in line.alt))
            assert len(line.sample) == 1
            maxref = max(maxref, len(line.ref))
            maxalt = max(maxalt, max((len(i) for i in line.alt)))

            if len(line.ref) != 1:
                print(str(line))
            if max((len(i) for i in line.alt)) != 1:
                print(str(line))
        assert maxref > 1
        assert maxalt > 1
