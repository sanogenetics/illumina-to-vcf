import itertools

import pytest

from illumina2vcf.illumina import IlluminaReader

from .conftest import IlluminaBuilder


class TestIllumina:
    @pytest.mark.parametrize("sep,sano", itertools.product(["\t", ","], [False, True]))
    def test_parse_header(self, sep: str, sano: bool) -> None:
        """
        GIVEN an illumina header
        """
        reader = IlluminaReader(sep)
        header = f"""[Header]
GSGT Version	2.0.4
Processing Date	2022-06-02 9:25 AM
Content		GSAMD-24v3-0-EA_20034606_A2.bpm
Num SNPs	730059
Total SNPs	730059
Num Samples	24
Total Samples	24"""
        if sano:
            header += "\nSANO"
        header = header.replace("\t", sep)
        """
        WHEN it is parsed
        """
        date, source = reader.parse_header(header.splitlines())

        """
        THEN it should have date and source as expected
        """
        # date must be in YYYYMMDD
        assert date == "20220602"

        # source must be GSAMD-24v3-0-EA_20034606_A2.bpm
        assert source == "GSAMD-24v3-0-EA_20034606_A2"

    @pytest.mark.parametrize("sano", (False, True))
    def test_generate_blocks(self, sano: bool) -> None:
        """
        GIVEN an illumina sorted file and reader
        """
        reader = IlluminaReader("\t")
        illumina = IlluminaBuilder().sano(sano).build_file()

        """
        WHEN it is parsed
        """
        _, _ = reader.parse_header(illumina)
        blocks = tuple(reader.generate_line_blocks(illumina))

        """
        THEN it should have blocks of lines
        """
        assert len(blocks) > 0
        for block in blocks:
            assert len(block) > 0

    @pytest.mark.parametrize("sano", (False, True))
    def test_block_sorting(self, sano: bool) -> None:
        """
        GIVEN an illumina unsorted file and reader
        """
        reader = IlluminaReader("\t")
        illumina = IlluminaBuilder().sano(sano).sorted(False).build_file()

        """
        WHEN it is parsed
        """
        _, _ = reader.parse_header(illumina)
        blocks = tuple(reader.generate_line_blocks(illumina))

        """
        THEN it should have blocks of lines
        """
        assert len(blocks) > 0
        for block in blocks:
            assert len(block) > 0

        # TODO better tests here, more types of block, etc
