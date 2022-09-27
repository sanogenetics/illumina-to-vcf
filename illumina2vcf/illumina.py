import csv
import logging
from typing import Dict, FrozenSet, Generator, Iterable, List, Tuple

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


class PositionError(Exception):
    pass


class HeaderError(Exception):
    pass


class DateError(Exception):
    pass


class DataError(Exception):
    pass


class IlluminaReader:
    delimiter: str
    blocklist_filename: str
    _blocklist: FrozenSet[str] = frozenset()

    def __init__(self, delimiter: str, blocklist_filename: str = ""):
        self.delimiter = delimiter
        self.blocklist_filename = blocklist_filename

    def __repr__(self) -> str:
        return f"IlluminaReader({repr(self.delimiter)},{repr(self.blocklist_filename)}"

    @property
    def blocklist(self) -> FrozenSet[str]:
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
            # handle files with a header line without a delimiter in them
            # e.g. legacy SANO files
            if self.delimiter in line:
                key, value = line.split(self.delimiter, 1)
            else:
                key = line
                value = ""

            file_header[key] = value
            # safety valve
            if len(file_header) >= 100:
                raise HeaderError("Header too long")

        if "Content" not in file_header:
            raise HeaderError("Unable to find Content field in header")

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
            raise DateError(f"Cannot parse Processing date '{date}' - not 3 components")
        if len(date_components[2]) == 4:
            # four digit year at end --- assume m/d/y
            date_components = [date_components[2], date_components[0], date_components[1]]
        elif len(date_components[0]) != 4:
            # four digit year at neither start nor end
            raise DateError(f"Cannot parse Processing date '{date}' - not 4 digit year")

        if int(date_components[1]) > 12:
            # not month in middle
            raise DateError(f"Cannot parse Processing date '{date}' - not 1-12 month")
        # ensure leading zeros
        date_components[1] = date_components[1].zfill(2)
        date_components[2] = date_components[2].zfill(2)

        date = "".join(date_components)
        return (date, source)

    def generate_line_blocks(self, input: Iterable[str]) -> Generator[List[Dict[str, str]], None, None]:
        block: List[dict] = []
        for i, row in enumerate(csv.DictReader(input, delimiter=self.delimiter)):
            if CHR not in row:
                raise DataError(f"{CHR} missing in row {i}")
            if POSITION not in row:
                raise DataError(f"{POSITION} missing in row {i}")
            if SNP_NAME not in row:
                raise DataError(f"{SNP_NAME} missing in row {i}")
            if SNP not in row:
                raise DataError(f"{SNP} missing in row {i}")

            # is there an existing block that has ended?
            if block and (block[0][CHR] != row[CHR] or block[0][POSITION] != row[POSITION]):
                # check we haven't gone backwards
                if block[0][CHR] == row[CHR] and int(block[0][POSITION]) > int(row[POSITION]):
                    raise PositionError(f"Position decrease on {row[CHR]} from {block[0][POSITION]} to {row[POSITION]}")
                # complete the previous block
                yield block
                block = []
            # ignore probes without chr or pos
            if row[CHR] == 0 or row[POSITION] == 0:
                logger.debug(f"Ignoring blocked {row[SNP_NAME]}")
                continue
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
