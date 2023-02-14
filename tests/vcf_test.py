import gzip
from typing import Dict, List, Tuple

from pytest import fixture, mark

from illumina2vcf import IlluminaReader, VCFMaker
from illumina2vcf.bpm.BPMReader import CSVManifestReader, ManifestFilter
from illumina2vcf.bpm.ReferenceGenome import ReferenceGenome

from .conftest import IlluminaBuilder


@fixture
def blocks() -> Tuple[List[Dict[str, str]], ...]:
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

    @mark.parametrize(
        "manifest_file",
        [
            open("tests/data/GSA-24v3-0_A2.trim.csv", "rt"),
            gzip.open("tests/data/GSA-24v3-0_A2.trim.csv.gz", "rt"),
        ],
        ids=["txt", "gz"],
    )
    def test_data_indel(self, blocks, genome_reader, manifest_file):
        """
        GIVEN an interable of blocks and indel data
        """
        indel_records = ManifestFilter(frozenset(), skip_snps=True).filtered_records(
            CSVManifestReader(manifest_file, genome_reader)
        )
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

        # pick out some specific ones that we know what to expect
        rsidlines = {
            "rs9651229": None,  # SNP
            "rs797044837": None,  # deletion
            "rs61750435": None,  # insertion
        }
        variant_info = {
            "rs9651229": "chr1-632287-C-T",
            "rs797044837": "chr1-1338000-CT-C",
            "rs61750435": "chr1-2406791-C-CT",
        }
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

            if line._id[0] in rsidlines.keys():
                rsid = line._id[0]
                assert variant_info[rsid] == "-".join([line.chrom, str(line.pos), line.ref, line.alt[0]])
                assert not rsidlines[rsid]
                rsidlines[rsid] = line
        assert maxref > 1
        assert maxalt > 1

        for key in rsidlines:
            assert rsidlines[key]
