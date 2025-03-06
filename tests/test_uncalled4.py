import copy

from pyfaidx import Fasta

from nanocompore.config import Config
from nanocompore.uncalled4 import Uncalled4

from tests.common import BASIC_CONFIG

def test_kmer_data_generator():
    config_yaml = copy.deepcopy(BASIC_CONFIG)
    config_yaml['data'] = {
        'cond1': {
            'sample1': {
                'bam': 'tests/fixtures/uncalled4_sample.bam'
            }
        },
        'cond2': { 
            'sample2': {
                'bam': 'tests/fixtures/uncalled4_sample.bam'
            }
        }
    }
    config_yaml['depleted_condition'] = 'cond1'
    config_yaml['kit'] = 'RNA004'
    config = Config(config_yaml)
    ref = 'ENST00000234590.10|ENSG00000074800.16|OTTHUMG00000001773.11|OTTHUMT00000497295.2|ENO1-201|ENO1|1781|protein_coding|'
    fasta_fh = Fasta('tests/fixtures/uncalled4_reference.fa')
    print(fasta_fh.keys())
    ref_seq = str(fasta_fh[ref])

    uncalled4 = Uncalled4(config, ref, ref_seq)
    kmers = list(uncalled4.kmer_data_generator())
    assert len(kmers) == 1715
    # The read has two signal segments aligned
    # to reference positions [6, 562) and [604, 1779).
    assert kmers[0].pos == 6
    assert kmers[0].kmer == 'TGTGGGTAC'
    # The last kmer of the first segment
    # should start at 553 (562 - 9).
    assert kmers[547].pos == 553
    assert kmers[547].kmer == 'AGTCCCGGC'
    # The first kmer of the second segment
    # starts at 604.
    assert kmers[548].pos == 604
    assert kmers[548].kmer == 'GGCCATGCA'
    # The last kmer starts at 1770.
    assert kmers[1714].pos == 1770
    assert kmers[1714].kmer == 'CCATGAGAA'

