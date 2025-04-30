import copy
import os
import shutil
import tempfile

from nanocompore.config import Config
from nanocompore.run import RunCmd
from tests.common import BASIC_CONFIG


def test_run():
    """
    A very simple integration test that just
    checks that the run command can finish
    successfully and produces the correct
    number of rows for the fixture bams.
    """
    config_yaml = copy.deepcopy(BASIC_CONFIG)
    output_path = tempfile.mkdtemp(prefix="test_")
    config_yaml['outpath'] = output_path
    config_yaml['min_coverage'] = 1
    try:
        config = Config(config_yaml)
        run_cmd = RunCmd(config)

        run_cmd()

        with open(os.path.join(output_path, "out_nanocompore_results.tsv")) as f:
            assert len(f.readlines()) == 3520
    finally:
        shutil.rmtree(output_path)

