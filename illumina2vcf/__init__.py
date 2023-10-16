import gzip
import logging
from typing import Iterable, TextIO, Union

from fsspec.core import OpenFile

from illumina2vcf.bpm.bpmreader import CSVManifestReader, ManifestFilter
from illumina2vcf.bpm.referencegenome import ReferenceGenome
from illumina2vcf.illumina import IlluminaReader
from illumina2vcf.vcf import VCFMaker

logger = logging.getLogger(__name__)


class Converter:
    reference: Union[str, OpenFile]
    reference_index: Union[str, OpenFile]
    blocklist_filename: str
    manifest_filename: str
    delimiter: str

    def __init__(
        self,
        reference: Union[str, OpenFile],
        reference_index: Union[str, OpenFile],
        blocklist_filename: str = "",
        manifest_filename: str = "",
        delimiter: str = ",",
        buildname="GRCh38",
    ) -> None:
        self.reference = reference
        self.reference_index = reference_index
        self.blocklist_filename = blocklist_filename
        self.manifest_filename = manifest_filename
        self.delimiter = delimiter
        self.buildname = buildname

    def convert(self, source: Iterable[str], destination: TextIO) -> None:
        reader = IlluminaReader(self.delimiter, self.blocklist_filename)

        genome_reader = ReferenceGenome(self.reference, self.reference_index)
        if self.manifest_filename:
            if self.manifest_filename.endswith(".gz"):
                logger.info(f"Reading gzip compressed manifest {self.manifest_filename}")
                manifest_reader = CSVManifestReader(gzip.open(self.manifest_filename, "rt"), genome_reader)
            else:
                logger.info(f"Reading uncompressed manifest {self.manifest_filename}")
                manifest_reader = CSVManifestReader(open(self.manifest_filename), genome_reader)
            indel_records = ManifestFilter(frozenset(), skip_snps=True).filtered_records(manifest_reader)
        else:
            indel_records = {}
        vcfgenerator = VCFMaker(genome_reader, indel_records)

        # read source header
        date, header_source = reader.parse_header(source)

        # read source blocks
        blocks = reader.generate_line_blocks(source)

        # generate vcf lines and qc info
        vcf_lines = list(vcfgenerator.generate_lines(blocks))

        # write header
        for line in vcfgenerator.generate_header(date, header_source, self.buildname):
            destination.write(str(line))
            destination.write("\n")

        # write vcf lines
        for line in vcf_lines:
            destination.write(str(line))
            destination.write("\n")
