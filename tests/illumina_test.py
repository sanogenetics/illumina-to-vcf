import pytest

from illumina2vcf.illumina import IlluminaReader


class TestIllumina:
    @pytest.mark.parametrize("sep", ["\t", ","])
    def test_parse_header(self, sep):
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
        date, source = reader.parse_header(header.splitlines())

        # date must be in YYYYMMDD
        assert date == "20220602"

        # source must be GSAMD-24v3-0-EA_20034606_A2.bpm
        assert source == "GSAMD-24v3-0-EA_20034606_A2"

    def test_generate_blocks(self):
        reader = IlluminaReader("\t")
        input = """
        """
        blocks = tuple(reader.generate_line_blocks(input.splitlines()))

        raise NotImplementedError()
