import pytest
from nanocompore.TxComp import *
from scipy.stats import combine_pvalues
import numpy as np
from unittest import mock
from nanocompore.SimReads import SimReads
from nanocompore.SampComp import SampComp


@pytest.fixture()
def fasta_file(tmp_path):
    fasta_file=tmp_path/"reference.fa"
    with open(fasta_file, 'w') as f:
        f.write('>Ref_001\nAATGGGAACGCTTAAACTTACCAAAGTGGTAGGGCATTATAACCCGATGAGACGTGTTACCCCTCAAAGGGCGGTACACCAACTTCAGACTGAGGTTCGATCTCGTGAAATGTCAGCAAGACCTCCACTCCGAAGCCATGAAGCGTTGCCGGGGATTAGCAAGCATCATGTTTACGAATTATATATGCAGACCGTCTTAACCTGGTCCCTAACTAATAATTTATAGTTCTTAGGACGCTCTTGATGATACCCGACGCCCGGCGATCTTATCATCTTTGGCCCCTTCTTCCGGATAGTGTAACGATCATAATTTTCTGCGGAACCTGAAGTTTGCTTTGAGCAAAACTGAGGAGGTTAGTTCTAATATCTGTGTCGCCAAAAACAGCACGATTCTCGACCCGGCGCCGCCACTCGCGACAGCCTTGGGTCCCAGATGCGAACTAATACTAACGGCCCCGACTGCAGAGAAATTTCGGCACACACTCATCTTTTTGCAGCCG\n')
        f.write('>Ref_002\nAATGAGGGATCGAAGCGA\n')
    return(str(fasta_file))

@pytest.fixture()
def nanopolishcomp_test_files(tmp_path, fasta_file):
    for rep in [1,2]:
        SimReads (
            fasta_fn=fasta_file,
            outpath=str(tmp_path),
            outprefix="control_rep"+str(rep),
            run_type = "RNA",
            intensity_mod_loc=0,
            intensity_mod_scale=0,
            dwell_mod_loc=0,
            dwell_mod_scale=0,
            mod_reads_freq=0,
            mod_bases_freq=0,
            mod_bases_type="A",
            pos_rand_seed=66,
            log_level="debug",
            overwrite=True)
    
    for rep in [1,2]:
        SimReads (
            fasta_fn=fasta_file,
            outpath=str(tmp_path),
            outprefix="mod_rep"+str(rep),
            run_type = "RNA",
            intensity_mod_loc=0.1,
            dwell_mod_loc=0.1,
            mod_reads_freq=0.5,
            mod_bases_freq=0.25,
            mod_bases_type="A",
            pos_rand_seed=66,
            log_level="debug",
            overwrite=True)

    fn_dict={"S1": {
                "R1": str(tmp_path / "control_rep1.tsv"), 
                "R2": str(tmp_path / "control_rep2.tsv")
                },
            "S2": {
                "R1": str(tmp_path / "mod_rep1.tsv"),
                "R2": str(tmp_path / "mod_rep2.tsv")
                }
            }
    return((fasta_file, fn_dict))

@pytest.mark.parametrize("method", ["GMM", "KS", "TT", "MW"])
@pytest.mark.parametrize("context", [2,3])
@pytest.mark.parametrize("context_weight", ["uniform", "harmonic"])
def test_1(nanopolishcomp_test_files, method, context, context_weight):
    fasta_file=nanopolishcomp_test_files[0]
    fn_dict=nanopolishcomp_test_files[1]
    s = SampComp(eventalign_fn_dict=fn_dict,
            outpath="tmp/",
            outprefix="nanocompore",
            fasta_fn=fasta_file,
            comparison_methods = method,
            logit = True,
            allow_warnings=True,
            sequence_context = context,
            sequence_context_weights = context_weight,
            min_coverage = 10,
            downsample_high_coverage = None,
            max_invalid_kmers_freq = 0.1,
            nthreads = 4,
            log_level = "debug",
            overwrite=True)
    
    db = s()
    
    db.save_report("tmp/report.txt")
