import csv
import logging
from typing import Dict, Generator, Iterable, List, Tuple

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

    def __repr__(self) -> str:
        return f"IlluminaReader({repr(self.delimiter)},{repr(self.blocklist_filename)}"

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

    def parse_header(self, header: Iterable[str]) -> Tuple[str, str]:
        file_header = {}
        # turn header lines into key/value
        for line in header:
            line = line.strip()
            if line == "[Header]":
                continue
            if line == "[Data]":
                break
            assert self.delimiter in line, f"'{self.delimiter}' not in '{line}'"
            key, value = line.split(self.delimiter, 1)
            file_header[key] = value
            # safety valve
            assert len(file_header) < 100

        source = file_header["Content"].split(".", 1)[0][1:].strip()  # remove trailing .pbm and leading sep char

        # need some validation on the dates here because there's a good chance they
        # will switch up the format on us at some point
        # this will currently work for month/day/year and year-month-day
        # (I guess also for month-day-year and year/month/day)
        # TODO use proper datetime parsing
        date = file_header["Processing Date"].strip().split(" ")[0]
        date_components = date.replace("/", "-").split("-")
        # not three pieces
        if len(date_components) != 3:
            raise DateError(f"Cannot parse Processing date '{date}' from line '{file_header[2]}'")

        if len(date_components[2]) == 4:
            # four digit year at end --- assume m/d/y
            date_components = [date_components[2], date_components[0], date_components[1]]
        elif len(date_components[0]) != 4:
            # four digit year at neither start nor end
            raise DateError(f"Cannot parse Processing date '{date}' from line '{file_header[2]}'")

        if int(date_components[1]) > 12:
            # not month in middle
            raise DateError(f"Cannot parse Processing date '{date}' from line '{file_header[2]}'")
        # ensure leading zeros
        date_components[1] = date_components[1].zfill(2)
        date_components[2] = date_components[2].zfill(2)

        date = "".join(date_components)
        return (date, source)

    def generate_line_blocks(self, input: Iterable[str]) -> Generator[List[Dict[str, str]], None, None]:
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
