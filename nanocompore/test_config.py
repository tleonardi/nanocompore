from nanocompore import Config


BASIC_VALID_CONFIG = {
    'data': {
        'kd': {
            'kd1': {
                'pod5': 'tests/fixtures/kd1.pod5',
                'bam': 'tests/fixtures/kd1.bam'
            },
            'kd2': {
                'pod5': 'tests/fixtures/kd2.pod5',
                'bam': 'tests/fixtures/kd2.bam'
            }
        },
        'wt': {
            'wt1': {
                'pod5': 'tests/fixtures/wt1.pod5',
                'bam': 'tests/fixtures/wt1.bam'
            },
            'wt2': {
                'pod5': 'tests/fixtures/wt2.pod5',
                'bam': 'tests/fixtures/wt2.bam'
            }
        }
    },
    'fasta': 'tests/fixtures/test_reference.fa',
}


class TestConfig:
    def test_valid_config(self):
       yaml = BASIC_VALID_CONFIG
       
       config = Config(yaml)
       
       assert config.getfasta_ref() == 'tests/fixtures/test_reference.fa'
       assert config.get_data() == yaml['data']
    
