from nanocompore.common import get_references_from_bams
from nanocompore.config import Config
from tests.common import BASIC_CONFIG


def test_get_references_from_bams():
    config = Config(BASIC_CONFIG)
    references = {t.ref_id for t in get_references_from_bams(config)}
    assert references == {'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|',
                          'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|'}

