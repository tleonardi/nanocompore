import os

import torch


cwd = os.getcwd()

BASIC_CONFIG = {
    'data': {
        'kd': {
            'kd1': {
                'bam': os.path.join(cwd, 'tests/fixtures/kd1.bam'),
            },
            'kd2': {
                'bam': os.path.join(cwd, 'tests/fixtures/kd2.bam'),
            }
        },
        'wt': {
            'wt1': {
                'bam': os.path.join(cwd, 'tests/fixtures/wt1.bam'),
            },
            'wt2': {
                'bam': os.path.join(cwd, 'tests/fixtures/wt2.bam'),
            }
        }
    },
    'fasta': os.path.join(cwd, 'tests/fixtures/test_reference.fa'),
    'depleted_condition': 'kd',
    'kit': 'RNA002',
    'resquiggler': 'uncalled4',
}

BASIC_EVENTALIGN_CONFIG = {
    'data': {
        'kd': {
            'kd1': {
                'db': os.path.join(cwd, 'tests/fixtures/kd1_eventalign.db'),
            },
            'kd2': {
                'db': os.path.join(cwd, 'tests/fixtures/kd2_eventalign.db'),
            }
        },
        'wt': {
            'wt1': {
                'db': os.path.join(cwd, 'tests/fixtures/wt1_eventalign.db'),
            },
            'wt2': {
                'db': os.path.join(cwd, 'tests/fixtures/wt2_eventalign.db'),
            }
        }
    },
    'fasta': os.path.join(cwd, 'tests/fixtures/test_reference.fa'),
    'depleted_condition': 'kd',
    'kit': 'RNA002',
    'resquiggler': 'eventalign',
}

class MockWorker:
    def log(self, level, msg):
        pass


def naneq(a, b):
    """
    Compare two torch tensors for equality, handling nans.
    """
    if a.shape != b.shape:
        return False
    nans = a.isnan() & b.isnan()
    return torch.eq(a[~nans], b[~nans]).all().item()

