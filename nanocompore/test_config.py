import pytest
import schema

import copy
import os

from nanocompore.Config import Config


cwd = os.getcwd()

BASIC_VALID_CONFIG = {
    'data': {
        'kd': {
            'kd1': {
                'pod5': os.path.join(cwd, 'tests/fixtures/kd1.pod5'),
                'bam': os.path.join(cwd, 'tests/fixtures/kd1.bam')
            },
            'kd2': {
                'pod5': os.path.join(cwd, 'tests/fixtures/kd2.pod5'),
                'bam': os.path.join(cwd, 'tests/fixtures/kd2.bam')
            }
        },
        'wt': {
            'wt1': {
                'pod5': os.path.join(cwd, 'tests/fixtures/wt1.pod5'),
                'bam': os.path.join(cwd, 'tests/fixtures/wt1.bam')
            },
            'wt2': {
                'pod5': os.path.join(cwd, 'tests/fixtures/wt2.pod5'),
                'bam': os.path.join(cwd, 'tests/fixtures/wt2.bam')
            }
        }
    },
    'fasta': os.path.join(cwd, 'tests/fixtures/test_reference.fa'),
}


class TestConfig:
    def test_valid_config(self):
       yaml = BASIC_VALID_CONFIG

       config = Config(yaml)

       assert config.get_fasta_ref() == os.path.join(cwd, 'tests/fixtures/test_reference.fa')
       assert config.get_data() == yaml['data']


    def test_too_many_conditions(self):
        yaml = copy.deepcopy(BASIC_VALID_CONFIG)
        yaml['data']['third_condition'] = yaml['data']['kd']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert 'Only two conditions allowed' == str(exception_info.value)


    def test_too_few_conditions(self):
        yaml = copy.deepcopy(BASIC_VALID_CONFIG)
        del yaml['data']['kd']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert 'Only two conditions allowed' == str(exception_info.value)


    def test_missing_fasta(self):
        yaml = copy.deepcopy(BASIC_VALID_CONFIG)
        del yaml['fasta']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "Missing key: 'fasta'" == str(exception_info.value)


    def test_invalid_fasta(self):
        yaml = copy.deepcopy(BASIC_VALID_CONFIG)
        yaml['fasta'] = yaml['data']['kd']['kd1']['pod5']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "Invalid fasta file" == str(exception_info.value)


    def test_too_few_threads(self):
        yaml = copy.deepcopy(BASIC_VALID_CONFIG)
        yaml['nthreads'] = 1

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "nthreads must be >= 3" == str(exception_info.value)


    def test_unsupported_comparison_method(self):
        yaml = copy.deepcopy(BASIC_VALID_CONFIG)
        yaml['comparison_methods'] = ['unsupported_method']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "did not validate 'unsupported_method'" in str(exception_info.value)


    def test_incorrect_sequence_context(self):
        yaml = copy.deepcopy(BASIC_VALID_CONFIG)
        yaml['sequence_context'] = 5

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "sequence_context must be >= 0 and <= 4" == str(exception_info.value)


    def test_incorrect_pvalue_threshold(self):
        yaml = copy.deepcopy(BASIC_VALID_CONFIG)
        yaml['pvalue_threshold'] = 5

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "pvalue_threshold must be >= 0 and <= 1" == str(exception_info.value)