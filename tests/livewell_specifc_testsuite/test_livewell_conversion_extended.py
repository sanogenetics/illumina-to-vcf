"""Expanded SNP/manifest, chromosome, genotype, and conversion invariants."""
from types import SimpleNamespace
import pytest
from illumina2vcf.illumina import IlluminaRow
from illumina2vcf.vcf import VCFMaker, ConverterError, Probe, QcStats

class Genome:
    def __init__(self, bases=None, names=('1','X','Y','M','chr2')):
        self.bases = bases or {('1',100):'A',('1',200):'C',('X',100):'A',('Y',100):'A',('M',100):'A',('chr2',100):'A'}
        self.reference_fasta = SimpleNamespace(faidx=SimpleNamespace(index={x:SimpleNamespace(rlen=100000) for x in names}))
    def get_reference_bases(self,chrom,start,end):
        return self.bases.get((chrom,start+1),'') if end==start+1 else ''

def row(sample='S1',snp='[A/G]',a='A',b='G',chrom='1',pos=100,probe='rs1',strand='+'):
    return IlluminaRow(sample,probe,chrom,pos,snp,a,b,strand)

def convert(rows,samples=None,genome=None,bpm=None):
    return VCFMaker(genome or Genome(), bpm)._line_block_to_vcf_line(rows, samples or ['S1'])

@pytest.mark.parametrize('input_chrom,output_chrom', [('1','1'),('chr1','1'),('X','X'),('chrX','X'),('XY','X'),('MT','M'),('chrM','M'),('2','chr2')])
def test_chromosome_alias_maps_to_reference(input_chrom,output_chrom):
    result=convert([row(chrom=input_chrom)],genome=Genome())
    assert result.chrom==output_chrom

@pytest.mark.parametrize('chrom', ['0','chr0','chrUn','GL000','23'])
def test_unavailable_chromosome_rejected(chrom):
    with pytest.raises(ConverterError):convert([row(chrom=chrom)])

@pytest.mark.parametrize('bases,snp,call,expected', [
    ({('1',100):'A'},'[A/G]',('G','A'),'0/1'),
    ({('1',100):'A'},'[A/G]',('G','G'),'1/1'),
    ({('1',100):'A'},'[A/G]',('A','A'),'0/0'),
    ({('1',100):'A'},'[A/G]',('-','-'),'./.'),
])
def test_all_basic_snp_genotypes(bases,snp,call,expected):
    v=convert([row(snp=snp,a=call[0],b=call[1])],genome=Genome(bases))
    assert v.sample[0]['GT']==expected

@pytest.mark.parametrize('alleles,expected', [(('0','1'),'0/1'),(('1','0'),'0/1'),(('2','1'),'1/2'),(('.','.'),'./.')])
def test_genotype_normalisation(alleles,expected):
    assert VCFMaker(Genome()).format_vcf_genotype(*alleles)==expected

def test_duplicate_probe_rows_do_not_change_genotype():
    # Convert a single probe row
    single = convert([row(a='A', b='G')])

    # Convert the same probe row repeated twice
    duplicate = convert([
        row(a='A', b='G'),
        row(a='A', b='G')
    ])

    # Duplicate identical rows must not change the genotype
    assert duplicate.sample[0]['GT'] == single.sample[0]['GT']

    # Both should produce the expected heterozygous genotype
    assert duplicate.sample[0]['GT'] == '0/1'

    # REF and ALT must also remain unchanged
    assert duplicate.ref == single.ref
    assert tuple(duplicate.alt) == tuple(single.alt)

def test_disagreeing_duplicate_probes_are_missing():
    v=convert([row(a='A',b='A'),row(a='G',b='G')]); assert v.sample[0]['GT']=='./.'

def test_absent_reference_base_rejected():
    with pytest.raises(ConverterError,match='Unable to get reference'):
        convert([row(pos=201)])

def test_mismatch_reference_base_rejected():
    with pytest.raises(ConverterError,match='Reference'):
        convert([row(snp='[C/G]',a='C',b='G')])

@pytest.mark.parametrize('invalid', ['A/G','[A/N]','[A/G','', '[a/g]'])
def test_bad_probe_notation_rejected(invalid):
    with pytest.raises(ConverterError,match='Unexpcted probes'):
        convert([row(snp=invalid)])

def test_reverse_strand_probe_alleles_are_complemented():
    # A/G on minus strand declares T/C on reference strand; use T as REF
    g=Genome({('1',100):'T'})
    v=convert([row(snp='[A/G]',a='T',b='C',strand='-')],genome=g)
    assert v.ref=='T' and tuple(v.alt)==('C',) and v.sample[0]['GT']=='0/1'

def test_multiallelic_alt_sorted_independently_of_input_order():
    # Two probes at same locus with different alternate bases; ref A
    v=convert([row(snp='[A/T]',a='A',b='T',probe='rsT'),row(snp='[A/G]',a='A',b='G',probe='rsG')])
    assert tuple(v.alt)==('G','T')

@pytest.mark.parametrize('assay', [-1,2,99])
def test_invalid_assay_type_is_rejected(assay):
    maker=VCFMaker(Genome())
    with pytest.raises(ValueError,match='assay_type'):
        maker._list_genotypes({'A','G'},('A','G'),Probe('rs1',assay,'A','G'))

@pytest.mark.parametrize('assay', [None,0,1])
def test_supported_assay_types_produce_genotype_candidates(assay):
    maker=VCFMaker(Genome())
    calls=maker._list_genotypes({'A','G'},('A','G'),Probe('rs1',assay,'A','G'))
    assert ('A','G') in calls

def test_no_call_gives_no_candidate_genotypes():
    maker=VCFMaker(Genome())
    assert maker._list_genotypes({'A','G'},('-','-'),Probe('rs1',None,'A','G'))==set()

def test_probe_calls_outside_allele_set_rejected():
    maker=VCFMaker(Genome())
    with pytest.raises(ValueError,match='no in allele list'):
        maker._list_genotypes({'A','G'},('A','T'),Probe('rs1',None,'A','G'))

def test_oversized_block_rejected_before_conversion():
    with pytest.raises(ConverterError,match='Oversized block'):
        convert([row() for _ in range(11)])

def test_filter_pass_if_at_least_one_sample_called():
    v = convert(
        [row('S1'), row('S2', a='-', b='-')],
        samples=['S1', 'S2']
    )

    fields = str(v).split('\t')

    assert fields[6] == 'PASS'
    assert v.sample[0]['GT'] == '0/1'
    assert v.sample[1]['GT'] == './.'

def test_filter_fail_if_all_samples_missing():
    v = convert([row(a='-', b='-')])

    fields = str(v).split('\t')

    assert fields[6] == 'FAIL'
    assert v.sample[0]['GT'] == './.'

def test_qc_counts_autosomal_y_and_x():
    qc=QcStats()
    for chrom,gt in [('1','0/1'),('1','./.'),('X','0/1'),('Y','0/0'),('M','0/0')]:
        line=SimpleNamespace(chrom=chrom,sample=[{'GT':gt}])
        qc._update_qc_stats(line)
    assert qc.vcf_lines_autosomal==2 and qc.called_lines_autosomal==1
    assert qc.heterozygous_lines_autosomal==1
    assert qc.called_lines_x==1 and qc.heterozygous_lines_x==1
    assert qc.vcf_lines_y==1 and qc.called_lines_y==1
    assert qc.call_rate('autosomal')==0.5

@pytest.mark.parametrize('chrom', ['autosomal','Y'])
def test_zero_denominator_qc_requires_explicit_handling(chrom):
    with pytest.raises(ZeroDivisionError): QcStats().call_rate(chrom)

