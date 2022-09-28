import logging
from typing import Union

from fsspec.core import OpenFiles

from .bpm.BPMReader import CSVManifestReader, ManifestFilter
from .bpm.ReferenceGenome import ReferenceGenome
from .illumina import IlluminaReader
from .vcf import VCFMaker

logger = logging.getLogger(__name__)


class Converter:
    reference: Union[str, OpenFiles]
    reference_index: Union[str, OpenFiles]
    blocklist_filename: str
    manifest_file: str
    delimiter: str

    def __init__(
        self,
        reference: Union[str, OpenFiles],
        reference_index: Union[str, OpenFiles],
        blocklist_filename: str = "",
        manifest_file: str = "",
        delimiter: str = ",",
        buildname="GRCh38",
    ) -> None:
        self.reference = reference
        self.reference_index = reference_index
        self.blocklist_filename = blocklist_filename
        self.manifest_file = manifest_file
        self.delimiter = delimiter
        self.buildname = buildname

    def convert(self, source, destination) -> None:
        reader = IlluminaReader(self.delimiter, self.blocklist_filename)
        genome_reader = ReferenceGenome(self.reference, self.reference_index)
        manifest_reader = CSVManifestReader(self.manifest_file, genome_reader, logger)
        indel_records = ManifestFilter(manifest_reader, "", True, logger).filtered_records()
        vcfgenerator = VCFMaker(genome_reader, indel_records)
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
