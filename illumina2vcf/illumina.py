import csv
import logging
from typing import Dict, Generator, List, Tuple

# header constants
SAMPLE_ID = "Sample ID"
SNP_NAME = "SNP Name"
CHR = "Chr"
POSITION = "Position"
SNP = "SNP"
ALLELE1 = "Allele1 - Plus"
ALLELE2 = "Allele2 - Plus"
STRAND = "Plus/Minus Strand"


logger = logging.getLogger(__name__)


class DateError(Exception):
    pass


class IlluminaReader:
    delimiter: str

    def __init__(self, delimiter):
        self.delimiter = delimiter

    def parse_header(self, source) -> Tuple[str, str]:
        file_header = []
        # read from source
        for line in source:
            line = line.strip()
            if line == "[Data]":
                break
            file_header.append(line)
            # safety valve
            assert len(file_header) < 100

        """[Header]
GSGT Version	2.0.4
Processing Date	2022-06-02 9:25 AM
Content		GSAMD-24v3-0-EA_20034606_A2.bpm
Num SNPs	730059
Total SNPs	730059
Num Samples	24
Total Samples	24"""

        source = file_header[3].split(self.delimiter)[-1].split(".")[0].lstrip()

        # need some validation on the dates here because there's a good chance they
        # will switch up the format on us at some point
        # this will currently work for month/day/year and year-month-day
        # (I guess also for month-day-year and year/month/day)
        # TODO use proper datetime parsing
        date = file_header[2].split(self.delimiter)[-1].lstrip().split(" ")[0]
        date_components = date.replace("/", "-").split("-")
        if len(date_components) != 3:
            raise DateError(f"Cannot parse Processing date {date} from line {file_header[2]}")
        if len(date_components[2]) == 4:
            date_components = [date_components[2], date_components[0], date_components[1]]
        elif len(date_components[0]) != 4:
            raise DateError(f"Cannot parse Processing date {date} from line {file_header[2]}")
        if int(date_components[1]) > 12:
            raise DateError(f"Cannot parse Processing date {date} from line {file_header[2]}")
        date_components[1] = date_components[1].zfill(2)
        date_components[2] = date_components[2].zfill(2)
        date = "".join(date_components)
        return (date, source)

    def generate_line_blocks(self, input) -> Generator[List[Dict[str, str]], None, None]:
        block = []
        for row in csv.DictReader(input, delimiter=self.delimiter):
            # is there an existing block that has ended?
            if block and (block[0][CHR] != row[CHR] or block[0][POSITION] != row[POSITION]):
                yield block
                block = []
            # ignore rows that are indels
            # ignore rows that are illumina duplicates
            if row[SNP] not in ("[D/I]", "[I/D]") and not "ilmndup" in row[SNP_NAME]:
                block.append(row)
        # if there is a final bloc, yield that too
        if block:
            yield block
