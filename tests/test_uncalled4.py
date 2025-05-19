import copy

import numpy as np
import pysam

from pyfaidx import Fasta

from nanocompore.config import Config
from nanocompore.uncalled4 import Uncalled4

from tests.common import BASIC_CONFIG

def test_get_data():
    config_yaml = copy.deepcopy(BASIC_CONFIG)
    # Note we use the same bam that contains
    # a single read so it doesn't matter
    # which read in the result we use for
    # the assertions.
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
    ref_seq = str(fasta_fh[ref])

    bams = {sample: pysam.AlignmentFile(sample_def['bam'], 'rb')
            for cond, cond_def in config.get_data().items()
            for sample, sample_def in cond_def.items()}
    uncalled4 = Uncalled4(ref, len(ref_seq), bams, config.get_kit())
    tensor, sample_ids, condition_ids = uncalled4.get_data()
    # 1781 positions on the transcript, 2 reads and 2 vars (intensity and dwell)
    assert tensor.shape == (1781, 2, 2)
    # The read has two signal segments aligned
    # to reference positions [6, 562) and [604, 1779).
    # Let's look at the alignment in the 3'->5'
    # direction:
    # The first reference position included in
    # the measurements is 1778. Uncalled4
    # assigns the measurement to the most
    # influential (central) position,
    # which for RNA004 is 5th of the 9-mer.
    # Hence, Uncalled4 will assign the first
    # value to reference position 1774. We want
    # to use the kmer start (in 5'->3') direction
    # so we would expect the last measurement
    # to be for position 1770.
    # That would be the first value in the uc
    # tag in the bam, because the uc tag
    # lists values in the order of sequencing
    # (3'->5' for RNA).
    assert tensor[1770, 0, 0] == 7992
    assert np.isnan(tensor[1771, 0, 0])

    # The last kmer in the first segment covers
    # positions [604, 612]. This means that
    # Uncalled4 will have the last value for
    # the segment assigned to 608, but we'll
    # consider this value to be for 604 because
    # we use the kmer start (in 5'->3' direction).
    assert tensor[604, 0, 0] == -5072
    # Uncalled4 has NaNs (-32768) for the remainng
    # 4 reference positions of the kmer.
    assert np.isnan(tensor[603, 0, 0])
    assert np.isnan(tensor[602, 0, 0])
    assert np.isnan(tensor[601, 0, 0])
    assert np.isnan(tensor[600, 0, 0])

    # Then we have 4 NaNs in the form of -32768
    # values for the first 4 positions of the
    # second segment (in 3'->5' direction), i.e.
    # positions [554, 557].
    assert np.isnan(tensor[557, 0, 0])
    assert np.isnan(tensor[556, 0, 0])
    assert np.isnan(tensor[555, 0, 0])
    assert np.isnan(tensor[554, 0, 0])

    assert tensor[553, 0, 0] == -10820

    # The most extreme 5' position reported
    # by Uncalled4 is 10, which for us will
    # be position 6.
    # This is the last value in the bam
    # (because they're listed in the sequencing
    # order, i.e 3'->5' for RNA).
    assert tensor[6, 0, 0] == 3802
    # Position 5 shouldn't have a value.
    assert np.isnan(tensor[5, 0, 0])


def test_get_data_no_reads():
    config_yaml = copy.deepcopy(BASIC_CONFIG)
    # Note we use the same bam that contains
    # no reads for both conditions
    config_yaml['data'] = {
        'cond1': {
            'sample1': {
                'bam': 'tests/fixtures/empty.bam'
            }
        },
        'cond2': { 
            'sample2': {
                'bam': 'tests/fixtures/empty.bam'
            }
        }
    }
    config_yaml['depleted_condition'] = 'cond1'
    config_yaml['kit'] = 'RNA004'
    config = Config(config_yaml)
    ref = 'ENST00000234590.10|ENSG00000074800.16|OTTHUMG00000001773.11|OTTHUMT00000497295.2|ENO1-201|ENO1|1781|protein_coding|'
    fasta_fh = Fasta('tests/fixtures/uncalled4_reference.fa')
    ref_seq = str(fasta_fh[ref])

    bams = {sample: pysam.AlignmentFile(sample_def['bam'], 'rb')
            for cond, cond_def in config.get_data().items()
            for sample, sample_def in cond_def.items()}
    uncalled4 = Uncalled4(ref, len(ref_seq), bams, config.get_kit())
    tensor, sample_ids, condition_ids = uncalled4.get_data()

    assert tensor.shape == (1781, 0, 2)
    assert len(sample_ids) == 0
    assert len(condition_ids) == 0

