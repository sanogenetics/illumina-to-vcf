import logging

from .illumina import IlluminaReader
from .vcf import VCFMaker

logger = logging.getLogger(__name__)


class Converter:
    reference_filename: str
    reference_index_filename: str
    blocklist_filename: str
    delimiter: str

    def __init__(
        self,
        reference_filename: str,
        reference_index_filename: str,
        blocklist_filename: str = "",
        delimiter: str = ",",
        buildname="GRCh38",
    ) -> None:
        self.reference_filename = reference_filename
        self.reference_index_filename = reference_index_filename
        self.blocklist_filename = blocklist_filename
        self.delimiter = delimiter
        self.buildname = buildname

    def convert(self, source, destination) -> None:
        reader = IlluminaReader(self.delimiter, self.blocklist_filename)
        vcfgenerator = VCFMaker(self.reference_filename, self.reference_index_filename)
        # read source header
        date, header_source = reader.parse_header(source)
        # write header
        for line in vcfgenerator.generate_header(date, header_source, self.buildname):
            destination.write(str(line))
            destination.write("\n")
        # read source blocks
        blocks = reader.generate_line_blocks(source)
        # write vcf lines
        for line in vcfgenerator.generate_lines(blocks):
            destination.write(str(line))
            destination.write("\n")
