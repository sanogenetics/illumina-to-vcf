import logging
from typing import Dict, FrozenSet, Iterable, Tuple

from .BPMRecord import BPMRecord, IndelSourceSequence
from .IlluminaBeadArrayFiles import RefStrand
from .ReferenceGenome import ReferenceGenome

logger = logging.getLogger(__name__)


class CSVManifestReader:
    _genome_reader: ReferenceGenome
    _required_columns = (
        "sourcestrand",
        "ilmnstrand",
        "name",
        "chr",
        "mapinfo",
        "refstrand",
        "sourceseq",
        "snp",
        "addressb_id",
        "allelea_probeseq",
    )
    _source_file: str

    def __init__(self, csv_source: Iterable[str], genome_reader: ReferenceGenome):
        """
        Initialize a manifest reader from a CSV file or other string iterable

            csv_filename    Path to the CSV manifest
            genome_reader   Reference genome reader
        """
        self._source = csv_source
        self._genome_reader = genome_reader

    def get_bpm_records(self):
        """
        Get BPM records from the reader

        Yields:
            BPMRecord: Next BPMRecord in the file
        Raises:
            Exception - Manifest is missing required column
        """
        in_data = False
        idx = -1
        for line in self._source:

            if line.startswith("IlmnID,"):
                in_data = True
                header = line.rstrip().lower().split(",")
                required_column2idx = {}
                for required_column in self._required_columns:
                    try:
                        required_column2idx[required_column] = header.index(required_column)
                    except:
                        raise Exception("Manifest is missing required column " + required_column)
                continue

            if line.startswith("[Controls]"):
                in_data = False
                continue

            if in_data:
                idx += 1
                bits = line.rstrip().split(",")
                (
                    source_strand,
                    ilmn_strand,
                    name,
                    chrom,
                    map_info,
                    ref_strand,
                    source_seq,
                    snp,
                    addressb_id,
                    probe_a,
                ) = [bits[required_column2idx[column]] for column in self._required_columns]

                if "D" in snp:
                    indel_source_sequence = IndelSourceSequence(source_seq)
                else:
                    indel_source_sequence = None

                assay_type = 0 if addressb_id == "" else 1
                try:
                    yield BPMRecord(
                        name,
                        "0",
                        probe_a,
                        chrom,
                        int(map_info),
                        snp,
                        RefStrand.from_string(ref_strand),
                        assay_type,
                        indel_source_sequence,
                        source_strand,
                        ilmn_strand,
                        self._genome_reader,
                        idx,
                    )
                except Exception as error:
                    logger.warning("Failed to process entry for record %s: %s", name, str(error))


class ManifestFilter:
    loci_to_filter: FrozenSet[str]
    skip_snps: bool

    def __init__(self, loci_to_filter: Iterable[str], skip_snps: bool = False):
        """
        Return a new ManifestFilter. Will skip records as specified in constructor
        as well as records with chromosome or mapping of zero.

            manifest_reader     The source of BPM records
            loci_to_filter      A set of record names to skip, may be None
            skip_snps           Skip SNPs
        """
        self.loci_to_filter = frozenset(loci_to_filter)
        self.skip_snps = skip_snps

    def filtered_records(self, manifest_reader: CSVManifestReader) -> Dict[Tuple[str, int], BPMRecord]:
        filtered_records = {}
        for record in manifest_reader.get_bpm_records():

            if record.chromosome == "0" or record.pos == 0:
                continue

            if self.loci_to_filter and record.name in self.loci_to_filter:
                continue

            if not record.is_indel() and self.skip_snps:
                continue

            chrom = record.chromosome
            # sex cleanup
            if chrom == "XX" or chrom == "XY":
                chrom = "X"
            # force chr prefix
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"

            # store it
            if (chrom, record.pos) not in filtered_records:
                filtered_records[(chrom, record.pos)] = []
            filtered_records[(chrom, record.pos)].append(record)

        logger.info(f"Got {len(filtered_records)} records")
        assert filtered_records
        return filtered_records
