"""Manifest-dependent allele mapping and indel safety tests with synthetic BPM records."""
from types import SimpleNamespace
import pytest
from illumina2vcf.vcf import VCFMaker, ConverterError, STRANDSWAP
from illumina2vcf.illumina import IlluminaRow

class Genome:
    def __init__(self):
        self.reference_fasta=SimpleNamespace(faidx=SimpleNamespace(index={'1':SimpleNamespace(rlen=100000)}))
    def get_reference_bases(self, chrom,start,end):
        return 'A' if (chrom,start,end)==('1',99,100) else 'C' if (chrom,start,end)==('1',98,99) else 'A'

def row(snp='[I/D]',a='I',b='D',name='indel',strand='+'):
    return IlluminaRow('S1',name,'1',100,snp,a,b,strand)

def manifest(snp='[I/D]',name='indel',assay_type=1,ref_strand=1,sequence='TT',deletion=False):
    rec=SimpleNamespace(name=name,snp=snp,assay_type=assay_type,ref_strand=ref_strand,
                        chromosome='1',pos=100,indel_source_sequence=sequence,is_deletion=deletion)
    rec.get_indel_source_sequences=lambda strand: ('',sequence,'')
    return rec

@pytest.mark.parametrize('indel,expected',[(False,('A','ATT',100)),(True,('CTT','C',99))])
def test_indel_reference_alt_and_position(indel,expected):
    maker=VCFMaker(Genome())
    assert maker.get_alleles_for_indel(manifest(deletion=indel))==expected

@pytest.mark.parametrize('is_del,genotype,expected',[
    (False,('I','I'),'1/1'),(False,('D','D'),'0/0'),(False,('I','D'),'0/1'),
    (True,('I','I'),'0/0'),(True,('D','D'),'1/1'),(True,('I','D'),'0/1'),
])
def test_indel_genotype_matches_reference_orientation(is_del,genotype,expected):
    assert VCFMaker(Genome()).convert_indel_genotype_to_vcf(genotype,is_del)==expected

@pytest.mark.xfail(strict=True, reason="Indel no-call is written as 0/0 (insertion) or 1/1 (deletion) with FILTER PASS")
@pytest.mark.parametrize('is_del', [False, True])
def test_indel_nocall_is_missing_through_converter(is_del):
    maker = VCFMaker(Genome(), {('1', 100): [manifest(deletion=is_del)]})
    v = maker._line_block_to_vcf_line([row(a='-', b='-')], ['S1'])
    assert v.sample[0]['GT'] == './.'


def test_indel_record_without_manifest_rejected():
    with pytest.raises(ConverterError,match='No BPM record for indel'):
        VCFMaker(Genome())._line_block_to_vcf_line([row()],['S1'])

def test_matching_indel_manifest_produces_valid_vcf():
    rec=manifest()
    line=VCFMaker(Genome(),{('1',100):[rec]})._line_block_to_vcf_line([row()],['S1'])
    assert (line.ref,tuple(line.alt),line.pos)==('A',('ATT',),100)
    assert line.sample[0]['GT']=='0/1'

def test_conflicting_indel_manifest_records_rejected():
    rec1=manifest(sequence='TT')
    rec2=manifest(name='indel2',sequence='GG')
    maker=VCFMaker(Genome(),{('1',100):[rec1,rec2]})
    with pytest.raises(ConverterError,match='Mismatched indel alleles'):
        maker._line_block_to_vcf_line([row()],['S1'])


def test_mixed_indel_and_snp_probes_cannot_be_combined():
    maker=VCFMaker(Genome())
    with pytest.raises(ConverterError,match='contains indels and SNPs'):
        maker._line_block_to_vcf_line([row(),row('[A/G]','A','G','snp')],['S1'])

def test_manifest_missing_probe_uses_row_fallback_warning(caplog):
    # BPM probe names do not match supplied row probe name.
    rec=manifest(snp='[A/G]',name='different',sequence=None)
    maker=VCFMaker(Genome(),{('1',100):[rec]})
    v=maker._line_block_to_vcf_line([row('[A/G]','A','G','not_in_bpm')],['S1'])
    assert v.sample[0]['GT']=='0/1'
    assert 'falling back to row-derived alleles' in caplog.text

def test_strand_complement_mapping_is_involutive():
    for allele in 'ACGTID':
        assert STRANDSWAP[STRANDSWAP[allele]]==allele

def test_indel_record_without_source_sequence_rejected():
    maker = VCFMaker(Genome(), {('1', 100): [manifest(sequence=None)]})
    with pytest.raises(ConverterError, match='No indel BPM record'):
        maker._line_block_to_vcf_line([row()], ['S1'])

def test_indel_row_at_snp_only_manifest_locus_raises_value_error():
    # SUSPECT: ValueError is not caught by generate_lines, so the whole conversion aborts
    maker = VCFMaker(Genome(), {('1', 100): [manifest(snp='[A/G]', sequence=None)]})
    with pytest.raises(ValueError, match='no in allele list'):
        maker._line_block_to_vcf_line([row()], ['S1'])