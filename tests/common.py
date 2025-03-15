import os


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

