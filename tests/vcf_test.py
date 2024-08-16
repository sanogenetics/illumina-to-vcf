import itertools
import re
import gzip
from typing import Dict, List, Tuple

from pytest import fixture, mark

from illumina2vcf import IlluminaReader, VCFMaker
from illumina2vcf.vcf import Probe, STRANDSWAP
from illumina2vcf.bpm.bpmreader import CSVManifestReader, ManifestFilter
from illumina2vcf.bpm.referencegenome import ReferenceGenome

from .conftest import IlluminaBuilder

@fixture
def blocks() -> Tuple[List[Dict[str, str]], ...]:
    reader = IlluminaReader("\t")
    illumina = IlluminaBuilder().build_file()
    _, _ = reader.parse_header(illumina)
    blocks = tuple(reader.generate_line_blocks(illumina))
    return blocks

@fixture
def genotypes() -> Dict[str, Dict[str, Tuple]]:
    genotypes = {
        'rs76584377': { # rs76584377-138_Inf2_C_T
            'Sample0': ('C', 'T'),
            'Sample1': ('C', 'T'),
            'Sample2': ('C', 'C'),
            'Sample3': ('C', 'C'),
            'Sample4': ('C', 'C'),
            'Sample5': ('C', 'C'),
            'Sample6': ('C', 'T'),
            'Sample7': ('C', 'T'),
            'Sample8': ('T', 'T'),
        },
        'rs76584377.1': { # rs76584377-138_Inf1_C_T
            'Sample0': ('C', 'T'),
            'Sample1': ('C', 'T'),
            'Sample2': ('C', 'C'),
            'Sample3': ('C', 'C'),
            'Sample4': ('C', 'C'),
            'Sample5': ('C', 'T'),
            'Sample6': ('C', 'T'),
            'Sample7': ('T', 'T'),
            'Sample8': ('T', 'T'),
        },
        'rs76584377.2': { # rs76584377-138_Inf1_C_G'
            'Sample0': ('C', 'C'),
            'Sample1': ('-', '-'),
            'Sample2': ('-', '-'),
            'Sample3': ('C', 'G'),
            'Sample4': ('G', 'G'),
            'Sample5': ('C', 'C'),
            'Sample6': ('-', '-'),
            'Sample7': ('G', 'G'),
            'Sample8': ('-', '-'),
        }
    }
    return genotypes

@fixture
def results() -> Dict[str, Tuple]:
    results = {'Sample0': ('C', 'T'), # C/T probes het, C/G probe has drop-out
            'Sample1': ('C', 'T'), # C/T probes het, C/G probe is no call
            'Sample2': ('-', '-'), # C/T probes hom ref, C/G probe is no call (could be CG)
            'Sample3': ('C', 'G'), # C/T probes hom ref, C/G probe is het (gen1 C/T probe is a dropout)
            'Sample4': ('-', '-'), # C/T probes hom ref, C/G probe is GG (gen1 probes conflict)
            'Sample5': ('-', '-'), # C/T probes conflict
            'Sample6': ('C', 'T'), # C/T probes are unambiguous without GC probe
            'Sample7': ('G', 'T'), # C/T both gen 1 probes have drop-outs
            'Sample8': ('T', 'T'), # don't need the G/C probe because gen2 T hom is unambiguous
            }
    return results

@fixture
def samples() -> Tuple:
    samples = ("Sample0", "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Sample7", "Sample8")
    return samples

@fixture
def genome_reader() -> ReferenceGenome:
    return ReferenceGenome(
        "tests/data/Homo_sapiens_assembly38.trim.fasta", "tests/data/Homo_sapiens_assembly38.trim.fasta.fai"
    )


class TestVCF:

    def test_header_only(self, genome_reader):
        """
        GIVEN a VCF generator
        """
        vcfgenerator = VCFMaker(genome_reader, {})

        """
        WHEN generating a header
        """
        lines = tuple(str(i) for i in vcfgenerator.generate_header("date", "source", "build"))
        """
        THEN it should produce valid header output
        """
        for line in lines:
            assert line
            assert line[0] == "#"

    def test_header(self, blocks, genome_reader):
        """
        GIVEN an interable of blocks
        """
        vcfgenerator = VCFMaker(genome_reader, {})

        """
        GIVEN a VCF generator
        """
        # generate vcf lines and qc info
        vcf_lines = [line for line in vcfgenerator.generate_lines(blocks)]

        """
        WHEN generating a header
        """
        lines = tuple(
            str(i) for i in vcfgenerator.generate_header("2022-03-25 9:20 AM", "GSAMD-24v3-0-EA_20034606_A2", "GRCh38")
        )
        """
        THEN it should produce valid header output
        """
        header_data = {}
        for line in lines:
            assert line
            assert line[:2] == "##"
            field, data = line[2:].split("=", maxsplit=1)
            if field in ("contig", "qc_stats"):
                data_dict = {}
                for stat in data[1:-1].split(","):
                    name, value = stat.split("=")
                    data_dict[name] = value
                if field == "contig":
                    try:
                        header_data[field].append(data_dict)
                    except KeyError:
                        header_data[field] = [
                            data_dict,
                        ]
                else:
                    header_data[field] = data_dict
            else:
                header_data[field] = data
        assert header_data["fileformat"] == "VCFv4.3"
        assert "filedate" in header_data  # this is not getting properly formatted
        assert header_data["source"] == '"GSAMD-24v3-0-EA_20034606_A2, Sano Genetics"'
        assert header_data["FILTER"] == '<ID=PASS,Description="All filters passed">'
        assert header_data["FORMAT"] == "<ID=GT,Number=1,Type=String,Description=Genotype>"
        assert len(header_data["contig"]) == 5
        for chrom_info in header_data["contig"]:
            assert int(chrom_info["length"]) == 3000000
            assert chrom_info["assembly"] == "GRCh38"
            assert chrom_info["ID"] in ("chr1", "chr2", "chr10", "chrX", "chrY")
        for stat in header_data["qc_stats"]:
            assert float(header_data["qc_stats"][stat]) <= 1

    def test_data(self, blocks, genome_reader):
        """
        GIVEN an interable of blocks
        """
        vcfgenerator = VCFMaker(genome_reader, {})
        """
        GIVEN a VCF generator
        """
        for _ in vcfgenerator.generate_header("date", "source", "build"):
            pass
        lines = vcfgenerator.generate_lines(blocks)
        """
        THEN it should produce valid output
        """
        maxref = 0
        maxalt = 0
        for line in lines:
            if line.comment_key or line.comment_raw:
                continue
            assert line.ref
            assert line.alt
            assert max(len(i) for i in line.alt)
            assert len(line.sample) == 1
            maxref = max(maxref, len(line.ref))
            maxalt = max(maxalt, max(len(i) for i in line.alt))
        assert maxref == 1
        assert maxalt == 1

    @mark.parametrize(
        "manifest_file",
        [
            open("tests/data/GSA-24v3-0_A2.trim.csv"),
            gzip.open("tests/data/GSA-24v3-0_A2.trim.csv.gz", "rt"),
        ],
        ids=["txt", "gz"],
    )
    def test_data_indel(self, blocks, genome_reader, manifest_file):
        """
        GIVEN an interable of blocks and indel data
        """
        bpm_records = ManifestFilter(frozenset(), skip_snps=False).filtered_records(
            CSVManifestReader(manifest_file, genome_reader)
        )
        vcfgenerator = VCFMaker(genome_reader, bpm_records)
        """
        GIVEN a VCF generator
        """
        for _ in vcfgenerator.generate_header("date", "source", "build"):
            pass
        lines = vcfgenerator.generate_lines(blocks)
        """
        THEN it should produce valid output
        """

        # pick out some specific ones that we know what to expect
        rsidlines = {
            "rs9651229": None,  # SNP
            "rs797044837": None,  # deletion
            "rs61750435": None,  # insertion
        }
        variant_info = {
            "rs9651229": "chr1-632287-C-T",
            "rs797044837": "chr1-1338000-CT-C",
            "rs61750435": "chr1-2406791-C-CT",
        }
        maxref = 0
        maxalt = 0
        for line in lines:
            if line.comment_key or line.comment_raw:
                continue
            assert line.ref
            assert line.alt
            assert max(len(i) for i in line.alt)
            assert len(line.sample) == 1
            maxref = max(maxref, len(line.ref))
            maxalt = max(maxalt, max(len(i) for i in line.alt))

            if line._id[0] in rsidlines.keys():
                rsid = line._id[0]
                assert variant_info[rsid] == "-".join([line.chrom, str(line.pos), line.ref, line.alt[0]])
                assert not rsidlines[rsid]
                rsidlines[rsid] = line
        assert maxref > 1
        assert maxalt > 1

        for key in rsidlines:
            assert rsidlines[key]


    def test_simplify_block(self, genotypes, results, samples, genome_reader) -> None:
        """
        GIVEN an illumina file and reader with the specified genotypes for each probe
        """

        reader = IlluminaReader("\t")
        illumina = IlluminaBuilder().samples(samples).genotypes(genotypes).build_file()

        """
        WHEN it is parsed and blocks are simplifed
        """
        manifest_file = open("tests/data/GSA-24v3-0_A2.trim.csv")

        bpm_records = ManifestFilter(frozenset(), skip_snps=False).filtered_records(
            CSVManifestReader(manifest_file, genome_reader)
        )
        vcfgenerator = VCFMaker(genome_reader, bpm_records)

        for _ in vcfgenerator.generate_header("date", "source", "build"):
            pass

        _, _ = reader.parse_header(illumina)
        blocks = tuple(reader.generate_line_blocks(illumina))

        """
        THEN it should have a single genotype per sample which matches the expected result
        """
        assert len(blocks) > 0

        for block in blocks:
            if block[0].chrom == "1" and block[0].pos == 1331911:
                locus_records =  bpm_records[('chr1', 1331911)]
                genotypes = vcfgenerator._simplify_block(block, locus_records)

                for sample in samples:
                    assert genotypes[0][sample] == results[sample]


    def test_combine_calls(self, genotypes, results, samples, genome_reader) -> None:
        """
        GIVEN an illumina file and reader with the specified genotypes for each probe
        """
        reader = IlluminaReader("\t")
        illumina = IlluminaBuilder().samples(samples).genotypes(genotypes).build_file()

        """
        WHEN it is parsed and probes are combined
        """
        manifest_file = open("tests/data/GSA-24v3-0_A2.trim.csv")

        bpm_records = ManifestFilter(frozenset(), skip_snps=False).filtered_records(
            CSVManifestReader(manifest_file, genome_reader)
        )
        vcfgenerator = VCFMaker(genome_reader, bpm_records)

        for _ in vcfgenerator.generate_header("date", "source", "build"):
            pass

        _, _ = reader.parse_header(illumina)
        blocks = tuple(reader.generate_line_blocks(illumina))

        """
        THEN combined genotypes should have the expected results
        """
        assert len(blocks) > 0

        for block in blocks:

            if block[0].chrom == "1" and block[0].pos == 1331911:
                genotypes = {}
                alleles = set()
                probes = {}
                for record in bpm_records[('chr1', 1331911)]:
                    match = re.match(r"\[([ATCGID]+)\/([ATCGID]+)\]", record.snp)
                    if not match:
                        msg = f"Unexpcted probes {row.chrom}:{row.pos} {probes}"
                        raise ConverterError(msg)

                    new_alleles = match.group(1, 2)

                    if record.ref_strand == "-":
                        new_alleles = tuple(STRANDSWAP[allele] for allele in new_alleles)
                    for allele in new_alleles:
                        alleles.add(allele)
                    probes[record.name] = Probe(
                        name=record.name,
                        assay_type=record.assay_type,
                        allele1=new_alleles[0],
                        allele2=new_alleles[1])

                for row in block:
                    sampleid = row.sample_id
                    new_calls = (row.allele1, row.allele2)
                    probe = probes[row.snp_name]

                    if sampleid in genotypes:
                        previous_genotypes = genotypes[sampleid]
                    else:
                        previous_genotypes = set()
                    genotypes[sampleid] = vcfgenerator._combine_calls(previous_genotypes, new_calls, alleles, probe)

                for sample in samples:
                    if len(genotypes[sample]) == 1:
                        genotype = list(genotypes[sample])[0]
                    else:
                        genotype = ('-', '-')
                    assert genotype == results[sample]
