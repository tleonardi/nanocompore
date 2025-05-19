import pytest
import schema

import copy
import os

from nanocompore.config import Config
from tests.common import BASIC_CONFIG, cwd


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
        yaml['fasta'] = yaml['data']['kd']['kd1']['bam']

        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        assert "Invalid fasta file" == str(exception_info.value)


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


    def test_db_required_for_eventalign(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        yaml['resquiggler'] = 'eventalign'

        # error when missing both
        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        # should be okay if eventalign_db is set
        for cond in yaml['data']:
            for sample in yaml['data'][cond]:
                yaml['data'][cond][sample]['db'] = 'some/path'
        config = Config(yaml)


    def test_db_required_for_remora(self):
        yaml = copy.deepcopy(BASIC_CONFIG)
        yaml['resquiggler'] = 'remora'

        # error when missing both
        with pytest.raises(schema.SchemaError) as exception_info:
            config = Config(yaml)

        # should be okay if eventalign_db is set
        for cond in yaml['data']:
            for sample in yaml['data'][cond]:
                yaml['data'][cond][sample]['db'] = 'some/path'
        config = Config(yaml)

