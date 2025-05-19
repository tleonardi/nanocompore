import copy
import multiprocessing
import os
import shutil
import tempfile

import numpy as np
import torch

from pyfaidx import Fasta

from nanocompore.common import INTENSITY_POS
from nanocompore.common import DWELL_POS
from nanocompore.common import MOTOR_DWELL_POS
from nanocompore.config import Config
from nanocompore.run import RunCmd
from nanocompore.run import Worker
from nanocompore.run import Uncalled4Worker
from nanocompore.run import GenericWorker
from nanocompore.transcript import Transcript
from tests.common import BASIC_CONFIG
from tests.common import BASIC_EVENTALIGN_CONFIG
from tests.common import naneq


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

    # This is a workaround for a very nasty issue
    # where running this integration test after
    # test_comparisons would cause it to crash with
    # a segmentation fault.
    # I've traced the problem and it appears that when
    # we run the GMM in the following tests:
    # - test_comparisons.test_gmm_test_split_single_component
    # - test_comparisons.test_gmm_test_split_two_components
    # then initialising any tensors in the run flow here
    # causes the segfault in the GOMP library.
    # Probably we're hitting some incompatibility
    # between the multiprocessing and the MPI used
    # by pytorch.
    # The workaround is to use spawn() to start the
    # worker processes, so that the torch memory
    # is not shared between the processes.
    # We then revert the change back to using
    # fork().
    multiprocessing.set_start_method("spawn", force=True)
    try:
        config = Config(config_yaml)
        run_cmd = RunCmd(config)

        run_cmd()

        with open(os.path.join(output_path, "out_nanocompore_results.tsv")) as f:
            assert len(f.readlines()) == 3520
    finally:
        shutil.rmtree(output_path)
        multiprocessing.set_start_method("fork", force=True)


def test_prepare_data_no_motor():
    config = Config(BASIC_CONFIG)

    worker = Worker(1, None, None, None, 0, 'cpu', config)

    data = np.array([
        # pos 0
        [[1., 1.],
         [1., 1.]],
        # pos 1
        [[2., 2.],
         [2., 2.]],
        # pos 2
        [[3., 3.],
         [3., 3.]],
    ])
    sample_ids = np.arange(3)
    condition_ids = np.array([0, 0, 1])

    tensor, samples, conditions, positions = worker._prepare_data(data, sample_ids, condition_ids)

    assert isinstance(tensor, torch.Tensor)
    assert tensor.shape == (3, 2, 2)
    assert tensor[0, 0, INTENSITY_POS] == 1.
    # The dwell is log transformed (and we add an epsilon
    # ot avoid division by zero when the value is 0).
    assert tensor[0, 0, DWELL_POS] == np.log10(1.0 + 1e-10)

    assert isinstance(samples, torch.Tensor)
    assert torch.equal(samples, torch.tensor([0, 1, 2]))

    assert isinstance(conditions, torch.Tensor)
    assert torch.equal(conditions, torch.tensor([0, 0, 1]))

    assert isinstance(positions, torch.Tensor)
    assert torch.equal(positions, torch.tensor([0, 1, 2]))


def test_prepare_data_with_motor():
    config_yaml = copy.deepcopy(BASIC_CONFIG)
    config_yaml['motor_dwell_offset'] = 2
    config = Config(config_yaml)

    worker = Worker(1, None, None, None, 0, 'cpu', config)

    data = np.array([
        # pos 0
        [[1., 1.],
         [1., 1.]],
        # pos 1
        [[2., 2.],
         [2., 2.]],
        # pos 2
        [[3., 3.],
         [3., 3.]],
    ])
    sample_ids = np.arange(3)
    condition_ids = np.array([0, 0, 1])

    tensor, samples, conditions, positions = worker._prepare_data(data, sample_ids, condition_ids)

    assert isinstance(tensor, torch.Tensor)
    assert tensor.shape == (3, 2, 3)
    assert tensor[0, 0, INTENSITY_POS] == 1.
    # The dwell is log transformed (and we add an epsilon
    # ot avoid division by zero when the value is 0).
    assert tensor[0, 0, DWELL_POS] == np.log10(1.0 + 1e-10)
    # The first position gets the motor from the 3rd,
    # and the remaining two positions have no motor dwell.
    assert tensor[0, 0, MOTOR_DWELL_POS] == np.log10(3.0 + 1e-10)
    assert tensor[1, 0, MOTOR_DWELL_POS].isnan()
    assert tensor[2, 0, MOTOR_DWELL_POS].isnan()

    assert isinstance(samples, torch.Tensor)
    assert torch.equal(samples, torch.tensor([0, 1, 2]))

    assert isinstance(conditions, torch.Tensor)
    assert torch.equal(conditions, torch.tensor([0, 0, 1]))

    assert isinstance(positions, torch.Tensor)
    assert torch.equal(positions, torch.tensor([0, 1, 2]))


def test_downsample():
    config = Config(BASIC_CONFIG)

    worker = Worker(1, None, None, None, 0, 'cpu', config)

    # covered positions per read are respectively:
    # 2, 1, 3, 0, 3, 2
    data = torch.tensor([
        # pos 0
        [[1., 1.],
         [np.nan, np.nan],
         [1., 1.],
         [np.nan, np.nan],
         [1., 1.],
         [1., 1.]],
        # pos 1
        [[2., 2.],
         [2., 2.],
         [2., 2.],
         [np.nan, np.nan],
         [1., 1.],
         [2., 2.]],
        # pos 2
        [[np.nan, np.nan],
         [np.nan, np.nan],
         [3., 3.],
         [np.nan, np.nan],
         [3., 3.],
         [np.nan, np.nan]],
    ])

    sample_ids = torch.arange(6)
    condition_ids = torch.tensor([0, 0, 0, 1, 1, 1])

    max_reads = 2
    r = worker._downsample(data, sample_ids, condition_ids, max_reads)
    data, samples, conditions = r

    assert data.shape == (3, 4, 2)
    assert torch.equal(samples, torch.tensor([0, 2, 4, 5]))
    assert torch.equal(conditions, torch.tensor([0, 0, 1, 1]))


def test_filter_processed_transcripts():
    config = Config(BASIC_CONFIG)

    worker = Worker(1, None, None, None, 0, 'cpu', config)

    # coverage per position is respectively:
    # cond 0
    # 2, 3, 1
    # cond 1
    # 2, 2, 1
    data = torch.tensor([
        # pos 0
        [[1., 1.],
         [np.nan, np.nan],
         [1., 1.],
         [np.nan, np.nan],
         [1., 1.],
         [1., 1.]],
        # pos 1
        [[2., 2.],
         [2., 2.],
         [2., 2.],
         [np.nan, np.nan],
         [1., 1.],
         [2., 2.]],
        # pos 2
        [[np.nan, np.nan],
         [np.nan, np.nan],
         [3., 3.],
         [np.nan, np.nan],
         [3., 3.],
         [np.nan, np.nan]],
    ])

    positions = torch.arange(3)
    conditions = torch.tensor([0, 0, 0, 1, 1, 1])

    min_cov = 2
    filtered_data, positions = worker._filter_low_cov_positions(data, positions, conditions, min_cov)

    assert filtered_data.shape == (2, 6, 2)
    assert naneq(filtered_data, data[:2])
    assert torch.equal(positions, torch.tensor([0, 1]))


def test_read_data_uncalled4():
    config = Config(BASIC_CONFIG)

    worker = Uncalled4Worker(1, None, None, None, 0, 'cpu', config)
    worker.setup()
    ref_id = 'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|'
    fasta_fh = Fasta(config.get_fasta_ref())
    ref_seq = str(fasta_fh[ref_id])

    transcript = Transcript(1, ref_id, ref_seq)

    data, samples, conditions = worker._read_data(transcript)

    assert data.shape == (2554, 12, 2)
    assert not np.isnan(data[301:303, :, :]).any()

    # Make sure that the reads are ordered.
    # This is necessary in order to guarantee
    # determinism of the GMM fitting (The order
    # of the data is important, because the k++
    # initialization selects random points to
    # initalize the GMM).
    expected_sample_ids = np.array([3, 3, 0, 1, 3, 2, 0, 1, 2, 1, 0, 2])
    assert np.array_equal(samples, expected_sample_ids)

    expected_condition_ids = np.array([1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1])
    assert np.array_equal(conditions, expected_condition_ids)


def test_read_data_eventalign():
    config = Config(BASIC_EVENTALIGN_CONFIG)

    worker = GenericWorker(1, None, None, None, 0, 'cpu', config)
    worker.setup()
    ref_id = 'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|'
    fasta_fh = Fasta(config.get_fasta_ref())
    ref_seq = str(fasta_fh[ref_id])

    transcript = Transcript(1, ref_id, ref_seq)

    data, samples, conditions = worker._read_data(transcript)
    print(data)
    print(samples, conditions)

    # In eventalign we have additional information
    # from the resquigglers that we can use to see
    # if a kmer is valid (we can check if there's
    # a mismatch between the model and reference kmer
    # or if the model kmer is NNNNN).
    # We then filter reads that have too many invalid
    # positions. Hence, here we have only 4 reads.
    assert data.shape == (2554, 4, 2)

    # Make sure that the samples are ordered.
    # This is necessary in order to guarantee
    # determinism of the GMM fitting (The order
    # of the data is important, because the k++
    # initialization selects random points to
    # initalize the GMM).
    # Note, the DB returns the reads in the same
    # order, so we can just order the samples
    # and we don't need to spend time ordering
    # by read or even read the read id.
    expected_sample_ids = np.array([0, 0, 1, 2])
    assert np.array_equal(samples, expected_sample_ids)

    expected_condition_ids = np.array([0, 0, 0, 1])
    assert np.array_equal(conditions, expected_condition_ids)


def test_get_transcripts_for_processing_insufficient_coverage():
    config = Config(BASIC_EVENTALIGN_CONFIG)
    run_cmd = RunCmd(config)

    refs = run_cmd._get_transcripts_for_processing()

    assert isinstance(refs, set)
    assert len(refs) == 0


def test_get_transcripts_for_processing_sufficient_coverage():
    config_yaml = copy.deepcopy(BASIC_CONFIG)
    config_yaml['min_coverage'] = 2
    config = Config(config_yaml)
    run_cmd = RunCmd(config)

    refs = run_cmd._get_transcripts_for_processing()
    refs = {ref.ref_id for ref in refs}

    assert isinstance(refs, set)
    assert len(refs) == 2
    assert refs == {'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|',
                    'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|'}


def test_get_transcripts_for_processing_selected_refs(tmpdir):
    p = tmpdir.join("selected_refs.txt")
    p.write("ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|\n")

    config_yaml = copy.deepcopy(BASIC_CONFIG)
    config_yaml['min_coverage'] = 2
    config_yaml['selected_refs'] = str(p)
    config = Config(config_yaml)
    run_cmd = RunCmd(config)

    refs = run_cmd._get_transcripts_for_processing()
    refs = {ref.ref_id for ref in refs}

    assert isinstance(refs, set)
    assert len(refs) == 1
    assert refs == {'ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|'}


def test_get_transcripts_for_processing_ignored_refs(tmpdir):
    p = tmpdir.join("ignored_refs.txt")
    p.write("ENST00000674681.1|ENSG00000075624.17|OTTHUMG00000023268|-|ACTB-219|ACTB|2554|protein_coding|\n")

    config_yaml = copy.deepcopy(BASIC_CONFIG)
    config_yaml['min_coverage'] = 2
    config_yaml['ignored_refs'] = str(p)
    config = Config(config_yaml)
    run_cmd = RunCmd(config)

    refs = run_cmd._get_transcripts_for_processing()
    refs = {ref.ref_id for ref in refs}

    assert isinstance(refs, set)
    assert len(refs) == 1
    assert refs == {'ENST00000642480.2|ENSG00000075624.17|OTTHUMG00000023268|OTTHUMT00000495153.1|ACTB-213|ACTB|2021|protein_coding|'}

