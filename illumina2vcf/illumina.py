import csv
import logging
from dataclasses import dataclass
from typing import FrozenSet, Generator, Iterable, List, Tuple

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


# with python 3.10 or higher, add slots=True here to save memory
@dataclass(frozen=True)
class IlluminaRow:
    sample_id: str
    snp_name: str
    chrom: str
    pos: int
    snp: str
    allele1: str
    allele2: str
    strand: str


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
        return f"IlluminaReader({self.delimiter!r},{self.blocklist_filename!r}"

    @property
    def blocklist(self) -> FrozenSet[str]:
        if not self.blocklist_filename:
            return frozenset()
        if not self._blocklist:
            blocklist = set()
            with open(self.blocklist_filename) as input_:
                for line in input_:
                    blocklist.add(line.strip())
            self._blocklist = frozenset(blocklist)
        return self._blocklist

    def parse_header(self, header: Iterable[str]) -> Tuple[str, str]:
        file_header = {}
        # turn header lines into key/value
        for line in header:
            line = line.strip()  # noqa: PLW2901
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
            if len(file_header) >= 100:  # noqa: PLR2004
                msg = "Header too long"
                raise HeaderError(msg)

        if "Content" not in file_header:
            msg = "Unable to find Content field in header"
            raise HeaderError(msg)

        source = file_header["Content"].split(".", 1)[0][1:].strip()  # remove trailing .pbm and leading sep char

        # need some validation on the dates here because there's a good chance they
        # will switch up the format on us at some point
        # this will currently work for month/day/year and year-month-day
        # (I guess also for month-day-year and year/month/day)
        # TODO use proper datetime parsing
        date = file_header["Processing Date"].strip().split(" ")[0]
        date_components = date.replace("/", "-").split("-")
        # not three pieces
        if len(date_components) != 3:  # noqa: PLR2004
            msg = f"Cannot parse Processing date '{date}' - not 3 components"
            raise DateError(msg)
        if len(date_components[2]) == 4:  # noqa: PLR2004
            # four digit year at end --- assume m/d/y
            date_components = [date_components[2], date_components[0], date_components[1]]
        elif len(date_components[0]) != 4:  # noqa: PLR2004
            # four digit year at neither start nor end
            msg = f"Cannot parse Processing date '{date}' - not 4 digit year"
            raise DateError(msg)

        if int(date_components[1]) > 12:  # noqa: PLR2004
            # not month in middle
            msg = f"Cannot parse Processing date '{date}' - not 1-12 month"
            raise DateError(msg)
        # ensure leading zeros
        date_components[1] = date_components[1].zfill(2)
        date_components[2] = date_components[2].zfill(2)

        date = "".join(date_components)
        return (date, source)

    def _sorted_lines(self, input_: Iterable[str]) -> Tuple[IlluminaRow, ...]:
        # this might get a bit big
        row_minis: List[IlluminaRow] = []
        for i, row in enumerate(csv.DictReader(input_, delimiter=self.delimiter)):
            if CHR not in row:
                msg = f"{CHR} missing in row {i}"
                raise DataError(msg)
            if POSITION not in row:
                msg = f"{POSITION} missing in row {i}"
                raise DataError(msg)
            if SNP_NAME not in row:
                msg = f"{SNP_NAME} missing in row {i}"
                raise DataError(msg)
            if SNP not in row:
                msg = f"{SNP} missing in row {i}"
                raise DataError(msg)

            row_mini = IlluminaRow(
                row[SAMPLE_ID],
                row[SNP_NAME],
                row[CHR],
                int(row[POSITION]),
                row[SNP],
                row[ALLELE1],
                row[ALLELE2],
                row[STRAND],
            )

            row_minis.append(row_mini)

        # this might take some time
        return tuple(sorted(row_minis, key=lambda r: (r.chrom, r.pos, r.sample_id)))

    def generate_line_blocks(self, input_: Iterable[str]) -> Generator[List[IlluminaRow], None, None]:
        block: List[IlluminaRow] = []
        for row in self._sorted_lines(input_):
            # is there an existing block that has ended?
            if block and (block[0].chrom != row.chrom or block[0].pos != row.pos):
                # check we haven't gone backwards
                if block[0].chrom == row.chrom and int(block[0].pos) > int(row.pos):
                    msg = f"Position decrease on {row.chrom} from {block[0].pos} to {row.pos}"
                    raise PositionError(msg)
                # complete the previous block
                yield block
                block = []
            # ignore probes without chr or pos
            if row.chrom == 0 or row.pos == 0:
                logger.debug(f"Ignoring blocked {row.snp_name}")
                continue
            # ignore blocked probes
            if row.snp_name in self.blocklist:
                logger.debug(f"Ignoring blocked {row.snp_name}")
                continue
            # ignore rows that are illumina duplicates
            if "ilmndup" in row.snp_name:
                logger.debug(f"Ignoring ilmndup {row.snp_name}")
                continue

            block.append(row)

        # if there is a final bloc, yield that too
        if block:
            yield block
