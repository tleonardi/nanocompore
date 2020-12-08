import pytest
from scipy.stats import combine_pvalues
import numpy as np
from unittest import mock
from loguru import logger

from nanocompore.SimReads import SimReads
from nanocompore.TxComp import *
from nanocompore.Whitelist import Whitelist
from nanocompore.common import *

# set logger lever
set_logger("debug")

@pytest.fixture(scope="module")
def fasta_file(tmp_path_factory):
    fasta_file = tmp_path_factory.mktemp("fasta") / "reference.fa"
    with open(str(fasta_file), 'w') as f:
        f.write('>Ref_001\n')
        f.write('A'*1000+'\n')
        f.write('>Ref_002\n')
        f.write('A'*50+'\n')
    return(str(fasta_file))

@pytest.fixture(scope="module")
def nanopolishcomp_test_files(tmp_path_factory, fasta_file):
    tmp_path=tmp_path_factory.mktemp("generated_data")
    SimReads(fasta_fn=fasta_file, outpath=str(tmp_path), outprefix="reads", overwrite=True)
    fn_dict={"S1": {
            "R1": str(tmp_path / "reads.tsv"),
            },
        "S2": {
            "R1": str(tmp_path / "reads.tsv"),
            }
        }
    return((fasta_file, fn_dict))

def test_def_len_filter(nanopolishcomp_test_files):
    fasta_file, fn_dict = nanopolishcomp_test_files
    whitelist = Whitelist(eventalign_fn_dict = fn_dict, fasta_fn = fasta_file, min_coverage = 100)
    assert len(whitelist.ref_reads) == 1

def test_no_len_filter(nanopolishcomp_test_files):
    fasta_file, fn_dict = nanopolishcomp_test_files
    whitelist = Whitelist(eventalign_fn_dict = fn_dict, fasta_fn = fasta_file, min_coverage = 100, min_ref_length=0)
    assert len(whitelist.ref_reads) == 2

def test_high_len_filter(nanopolishcomp_test_files):
    fasta_file, fn_dict = nanopolishcomp_test_files
    whitelist = Whitelist(eventalign_fn_dict = fn_dict, fasta_fn = fasta_file, min_coverage = 100, min_ref_length=1000)
    assert len(whitelist.ref_reads) == 0

def test_high_cov_filter(nanopolishcomp_test_files):
    fasta_file, fn_dict = nanopolishcomp_test_files
    whitelist = Whitelist(eventalign_fn_dict = fn_dict, fasta_fn = fasta_file, min_coverage = 101)
    assert len(whitelist.ref_reads) == 0
