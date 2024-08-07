import itertools
from typing import Tuple

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
        header = """[Header]
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

    @pytest.mark.parametrize(
        "sano,sorted,data_header,samples",
        itertools.product(
            (False, True),
            (True, False),
            (True, False),
            (("sample1",), ("sample1", "sample2")),
        ),
        ids=[
            "-".join(i)
            for i in itertools.product(
                ["GSA", "Sano"],
                ["sorted", "unsorted"],
                ["fullheader", "miniheader"],
                ["1sample", "2sample"],
            )
        ],
    )
    def test_generate_blocks(self, sano: bool, sorted: bool, data_header: bool, samples: Tuple[str, ...]) -> None:
        """
        GIVEN an illumina file and reader
        """
        reader = IlluminaReader("\t")
        illumina = IlluminaBuilder().samples(samples).sano(sano).sorted(sorted).data_header(data_header).build_file()

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
            # each block should be fully divisible by the number of samples
            assert len(block) % len(samples) == 0
            # all lines in a block should have the same chromosome and position
            assert len(set(i.chrom for i in block)) == 1
            assert len(set(i.pos for i in block)) == 1
            # each line should have a unique sample and ilmn ids
            assert len(set((i.sample_id, i.ilmn_id) for i in block)) == len(block)
