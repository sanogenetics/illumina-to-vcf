import pytest

from illumina2vcf.illumina import IlluminaReader


class TestIllumina:
    @pytest.mark.parametrize("sep", ["\t", ","])
    def test_parse_header(self, sep):
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

    def test_generate_blocks(self, fixture_illumina_lines):
        """
        GIVEN an illumina file and reader
        """
        reader = IlluminaReader("\t")

        """
        WHEN it is parsed
        """
        date, header_source = reader.parse_header(fixture_illumina_lines)
        blocks = tuple(reader.generate_line_blocks(fixture_illumina_lines))

        """
        THEN it should have blocks of lines
        """
        assert len(blocks) > 0
        for block in blocks:
            assert len(block) > 0

        # TODO better tests here, more types of block, etc
