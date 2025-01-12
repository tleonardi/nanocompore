import pytest
import schema

import copy
import os

from nanocompore.config import Config


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
    'preprocessing_db': 'tests/fixtures/kmer_data.sqlite'
}


class TestConfig:
    def test_valid_config(self):
       yaml = BASIC_CONFIG

       config = Config(yaml)

       assert config.get_fasta_ref() == os.path.join(cwd, 'tests/fixtures/test_reference.fa')
       assert config.get_data() == yaml['data']


    def test_too_many_conditions(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        yaml['data']['third_condition'] = yaml['data']['kd']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert 'Only two conditions allowed' == str(exception_info.value)


    def test_too_few_conditions(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        del yaml['data']['kd']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert 'Only two conditions allowed' == str(exception_info.value)


    def test_missing_fasta(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        del yaml['fasta']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "Missing key: 'fasta'" == str(exception_info.value)


    def test_invalid_fasta(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        yaml['fasta'] = yaml['data']['kd']['kd1']['pod5']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "Invalid fasta file" == str(exception_info.value)


    def test_too_few_threads(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        yaml['nthreads'] = 1

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "nthreads must be >= 2" == str(exception_info.value)


    def test_unsupported_comparison_method(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        yaml['comparison_methods'] = ['unsupported_method']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "did not validate 'unsupported_method'" in str(exception_info.value)


    def test_comparison_method_is_case_sensitive(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        yaml['comparison_methods'] = ['gmm']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "did not validate 'gmm'" in str(exception_info.value)


    def test_incorrect_sequence_context(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        yaml['sequence_context'] = 5

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "sequence_context must be >= 0 and <= 4" == str(exception_info.value)


    def test_incorrect_pvalue_threshold(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        yaml['pvalue_threshold'] = 5

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "pvalue_threshold must be >= 0 and <= 1" == str(exception_info.value)


    def test_unsupported_correction_method(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        yaml['correction_method'] = 'bonferroni'

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "does not match 'bonferroni'" in str(exception_info.value)


    def test_missing_kit(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        del yaml['kit']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "Missing key: 'kit'" in str(exception_info.value)


    def test_unsupported_kit(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        yaml['kit'] = 'RNA500'

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "did not validate 'RNA500'" in str(exception_info.value)


    def test_missing_resquiggler(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        del yaml['resquiggler']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "Missing key: 'resquiggler'" in str(exception_info.value)


    def test_missing_depleted_condition(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        del yaml['depleted_condition']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "Missing key: 'depleted_condition'" in str(exception_info.value)


    def test_nonexisting_depleted_condition(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        yaml['depleted_condition'] = 'does_not_exist'

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)


    def test_bam_required_for_uncalled4(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        del yaml['data']['kd']['kd1']['bam']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)


    def test_eventalign_tsv_or_db_required_for_eventalign(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        yaml['resquiggler'] = 'eventalign'

        # should be okay if eventalign_tsv is set
        config = Config(yaml)

        file = yaml['data']['kd']['kd1']['eventalign_tsv']

        # error when missing both
        del yaml['data']['kd']['kd1']['eventalign_tsv']
        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        # should be okay if eventalign_db is set
        yaml['data']['kd']['kd1']['eventalign_db'] = file
        config = Config(yaml)

