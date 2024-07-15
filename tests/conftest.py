import csv
import random
from dataclasses import dataclass
from datetime import datetime
from io import StringIO
from typing import Any, Dict, Generator, Iterable, List, Tuple

from illumina2vcf.bpm.bpmrecord import COMPLEMENT_MAP


@dataclass
class Probe:
    ilmn_id: str
    name: str
    chrom: str
    pos: int
    a: str
    b: str
    strand: str


class IlluminaBuilder:
    _sano: bool = False
    _sorted: bool = True
    _rng: random.Random
    _data_header_full: bool = True
    _samples: Tuple[str, ...]
    _genotypes: Dict[str, Dict[str, Tuple[str, str]]]

    def __init__(self):
        self._rng = random.Random(42)
        self._samples = ("sample1",)
        self._genotypes = ()
        pass

    def sano(self, sano: bool) -> "IlluminaBuilder":
        self._sano = sano
        return self

    def sorted(self, sorted: bool) -> "IlluminaBuilder":
        self._sorted = sorted
        return self

    def data_header(self, data_header_full: bool) -> "IlluminaBuilder":
        self._data_header_full = data_header_full
        return self

    def samples(self, samples: Iterable[str]) -> "IlluminaBuilder":
        self._samples = tuple(samples)
        return self

    def genotypes(self, genotypes: Iterable[str]) -> "IlluminaBuilder":
        self._genotypes = genotypes
        return self

    def _generate_probes(self) -> Generator[Probe, None, None]:
        with open("tests/data/GSA-24v3-0_A2.trim.csv") as infile:
            # read through the header
            for line in infile:
                if line.strip() == "[Assay]":
                    break
            # switch to CSV
            reader = csv.DictReader(infile)
            for row in reader:
                yield Probe(
                    ilmn_id = row["IlmnID"],
                    name=row["Name"],
                    chrom=row["Chr"],
                    pos=int(row["MapInfo"]),
                    a=row["SNP"][1],  # [X/.]
                    b=row["SNP"][-2],  # [./X]
                    strand=row["RefStrand"],
                )

    def _generate_sorted_probes(self) -> List[Probe]:
        def probekey(probe: Probe) -> Tuple[str, int, int]:
            if probe.chrom.isnumeric():
                return ("", int(probe.chrom), probe.pos)
            else:
                return (probe.chrom, 0, probe.pos)

        return sorted(self._generate_probes(), key=probekey)

    def _generate_unsorted_probes(self) -> List[Probe]:
        # convert positions to string before sorting so that they will be out of
        # order (this assumes that there are positions with different numbers of digits)
        def str_probekey(probe: Probe) -> Tuple[str, int, str]:
            if probe.chrom.isnumeric():
                return ("", int(probe.chrom), str(probe.pos))
            else:
                return (probe.chrom, 0, str(probe.pos))

        return sorted(self._generate_probes(), key=str_probekey)

    def _generate_header_lines(self, num_snps=730059, num_samples=24):
        yield "[Header]"
        yield "GSGT Version	2.0.4"
        # e.g. "Processing Date	2022-06-02 9:25 AM"
        yield f"Processing Date	{datetime.now().strftime('%Y-%m-%d %I:%M %p')}"
        yield "Content		GSAMD-24v3-0-EA_20034606_A2.bpm"
        yield f"Num SNPs	{num_snps}"
        yield f"Total SNPs	{num_snps}"
        yield f"Num Samples	{num_samples}"
        yield f"Total Samples	{num_samples}"
        if self._sano:
            yield "SANO"

    @staticmethod
    def _map_to_line(header: Iterable[str], data: Dict[str, Any]) -> str:
        line_items = []
        for heading in header:
            line_items.append(str(data.get(heading, ".")))
        return "\t".join(line_items)

    def _generate_data_header(self) -> Tuple[str, ...]:
        if self._data_header_full:
            return tuple(
                (
                    "Sample ID",
                    "IlmnID",
                    "RsID",
                    "GC Score",
                    "SNP Name",
                    "SNP Index",
                    "Sample Index",
                    "Sample Name",
                    "Sample Group",
                    "SNP Aux",
                    "Chr",
                    "Position",
                    "GT Score",
                    "Cluster Sep",
                    "SNP",
                    "ILMN Strand",
                    "Customer Strand",
                    "Top Genomic Sequence",
                    "Plus/Minus Strand",
                    "Allele1 - Plus",
                    "Allele2 - Plus",
                    "Allele1 - Forward",
                    "Allele2 - Forward",
                    "Theta",
                    "R",
                    "X",
                    "Y",
                    "X Raw",
                    "Y Raw",
                    "B Allele Freq",
                    "Log R Ratio",
                    "CNV Value",
                    "CNV Confidence",
                )
            )
        else:
            return tuple(
                (
                    "Sample ID",
                    "IlmnID",
                    "SNP Name",
                    "Sample Name",
                    "Chr",
                    "Position",
                    "SNP",
                    "Plus/Minus Strand",
                    "Allele1 - Plus",
                    "Allele2 - Plus",
                )
            )

    def _generate_data_lines(self, probes, genotypes = None) -> Generator[str, None, None]:
        header = self._generate_data_header()
        yield "[Data]"
        yield "\t".join(header)
        for i, sample in enumerate(self._samples, 1):
            # TODO mark some as chr/pos 0/0
            for probe in probes:
                data = {}
                data["Sample ID"] = sample
                data["Sample Index"] = i
                data["IlmnID"] = probe.ilmn_id
                data["SNP Name"] = probe.name
                data["Chr"] = probe.chrom
                data["Position"] = probe.pos
                data["SNP"] = f"[{probe.a}/{probe.b}]"
                data["Plus/Minus Strand"] = probe.strand
                alleles = (
                    [COMPLEMENT_MAP[probe.a], COMPLEMENT_MAP[probe.b]] if probe.strand == "-" else [probe.a, probe.b]
                )
                if probe.ilmn_id in self._genotypes:
                    allele1 = self._genotypes[probe.ilmn_id][sample][0]
                    assert allele1 in alleles or allele1 == '-'
                    data["Allele1 - Plus"] = allele1
                    allele2 = self._genotypes[probe.ilmn_id][sample][1]
                    assert allele2 in alleles or allele2 == '-'
                    data["Allele2 - Plus"] = allele2
                else:
                    data["Allele1 - Plus"] = self._rng.choice(alleles)
                    data["Allele2 - Plus"] = self._rng.choice(alleles)
                yield self._map_to_line(header, data)

    def _generate_lines(self) -> Generator[str, None, None]:
        # first generate the probes
        # make sure they are sorted by chromosome and position
        if self._sorted:
            probes = self._generate_sorted_probes()
        else:
            probes = self._generate_unsorted_probes()
        # now we know how many probes and how many samples we have can generate header
        for line in self._generate_header_lines(num_snps=len(probes), num_samples=len(self._samples)):
            yield line

        # now we can generate data lines for each sample
        for line in self._generate_data_lines(probes):
            yield line

    def build_file(self) -> StringIO:
        return StringIO("\n".join(self._generate_lines()))
