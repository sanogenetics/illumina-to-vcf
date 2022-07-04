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
    blocklist_filename: str
    _blocklist: frozenset[str] = frozenset()

    def __init__(self, delimiter: str, blocklist_filename: str = ""):
        self.delimiter = delimiter
        self.blocklist_filename = blocklist_filename

    @property
    def blocklist(self) -> frozenset[str]:
        if not self.blocklist_filename:
            return frozenset()
        if not self._blocklist:
            blocklist = set()
            with open(self.blocklist_filename, "rt") as input:
                for line in input:
                    blocklist.add(line.strip())
            self._blocklist = frozenset(blocklist)
        return self._blocklist

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
        dateline = file_header[2]
        _,date = dateline.split(self.delimiter, 1)
        date = date.strip()
        date,_ = date.split(" ",1)
        date_components = date.replace("/", "-").split("-")
        if len(date_components) != 3:
            raise DateError(f"Cannot parse Processing date '{date}' from line '{file_header[2]}' - not 3 components")
        if len(date_components[2]) == 4:
            date_components = [date_components[2], date_components[0], date_components[1]]
        elif len(date_components[0]) != 4:
            raise DateError(f"Cannot parse Processing date '{date}' from line '{file_header[2]}' - not 4 digit year")
        if int(date_components[1]) > 12:
            raise DateError(f"Cannot parse Processing date '{date}' from line '{file_header[2]}' - not 1-12 month")
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
            # ignore blocked probes
            if row[SNP_NAME] in self.blocklist:
                logger.debug(f"Ignoring blocked {row[SNP_NAME]}")
                continue
            # ignore rows that are illumina duplicates
            if "ilmndup" in row[SNP_NAME]:
                logger.debug(f"Ignoring ilmndup {row[SNP_NAME]}")
                continue
            # ignore rows that are indels
            if row[SNP] in ("[D/I]", "[I/D]"):
                logger.debug(f"Ignoring indel {row[SNP_NAME]}")
                continue

            block.append(row)

        # if there is a final bloc, yield that too
        if block:
            yield block
