import os


cwd = os.getcwd()

BASIC_CONFIG = {
    'data': {
        'kd': {
            'kd1': {
                'pod5': os.path.join(cwd, 'tests/fixtures/kd1.pod5'),
                'bam': os.path.join(cwd, 'tests/fixtures/kd1.bam'),
                'eventalign_tsv': os.path.join(cwd, 'tests/fixtures/kd1_eventalign.tsv'),
            },
            'kd2': {
                'pod5': os.path.join(cwd, 'tests/fixtures/kd2.pod5'),
                'bam': os.path.join(cwd, 'tests/fixtures/kd2.bam'),
                'eventalign_tsv': os.path.join(cwd, 'tests/fixtures/kd2_eventalign.tsv'),
            }
        },
        'wt': {
            'wt1': {
                'pod5': os.path.join(cwd, 'tests/fixtures/wt1.pod5'),
                'bam': os.path.join(cwd, 'tests/fixtures/wt1.bam'),
'eventalign_tsv': os.path.join(cwd, 'tests/fixtures/wt1_eventalign.tsv'),
            },
            'wt2': {
                'pod5': os.path.join(cwd, 'tests/fixtures/wt2.pod5'),
                'bam': os.path.join(cwd, 'tests/fixtures/wt2.bam'),
                'eventalign_tsv': os.path.join(cwd, 'tests/fixtures/wt2_eventalign.tsv'),
            }
        }
    },
    'fasta': os.path.join(cwd, 'tests/fixtures/test_reference.fa'),
    'depleted_condition': 'kd',
    'kit': 'RNA004',
    'resquiggler': 'uncalled4',
    'preprocessing_db': os.path.join(cwd, 'tests/fixtures/kmer_data.sqlite')
}
