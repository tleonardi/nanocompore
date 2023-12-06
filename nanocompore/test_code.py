#!/usr/bin/env python 
from crypt import methods
import sqlite3, pytest
import nanocompore.Data_DB_manager as db_mng
import nanocompore.Eventalign_DB as db
import nanocompore.SampComp_SQLDB as results_db
import Nanocompore_rewrite.nanocompore.Transcript as Transcript
import nanocompore.Sample as Sample
import nanocompore.TxComp as TxComp
import nanocompore.SampComp as SampComp
import nanocompore.SampCompResultsmanager as SampCompResultsmanager
from nanocompore.common import *
import pandas as pd
import numpy as np

from loguru import logger

import os

sample2condition = {'KO1':"KO", 'KO2':"KO", 'KO3':"KO", 
                    'WT1':"WT", 'WT2':"WT", 'WT3':"WT"}
#test dataset
test_samples = {'KO1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/sqlite_KO1_eventalign_collapse.db',
                'WT1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/WT1_eventalign_collapse.db'}
#list of transcripts for testing purposes
tx = ['YCR010C_mRNA', 'YNL208W_mRNA', 'YNL037C_mRNA',
      'YDR099W_mRNA', 'YMR303C_mRNA', 'YNL155W_mRNA',
      'YGL009C_mRNA', 'YEL060C_mRNA', 'YBR067C_mRNA',
      'YPR036W-A_mRNA', 'YCR005C_mRNA']

test_tx_name = tx[0]
tx_length = 852
max_coverage = float('inf')
min_coverage = 10
pos_0 = 1
test_methods = ['KS', 'GMM']

bed_fn = ''
fasta_fn = "/nfs/research/birney/users/logan/refs/yeast/demo_yeast.fa"

outpath = '/hps/nobackup/birney/users/logan/nanocompore_rewrite/'
prefix = ''

if prefix:
    prefix = f"{prefix}_"

subcommand = 'SampComp'

log_fn = os.path.join(outpath, f"{prefix}{subcommand}.log")
log_level='info'

set_logger(log_level, log_fn=log_fn)

def one_sample_per_condition():
    test_reference_samples = {'KO1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/sqlite_KO1_eventalign_collapse.db'}

    test_test_samples = {'WT1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/WT1_eventalign_collapse.db'}

    return test_reference_samples, test_test_samples

def two_sample_per_condition():
    test_reference_samples = {'KO1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/sqlite_KO1_eventalign_collapse.db',
                              'KO2':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/sqlite_KO1_eventalign_collapse2.db'}

    test_test_samples = {'WT1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/WT1_eventalign_collapse.db',
                         'WT2':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/WT1_eventalign_collapse2.db'}

    return test_reference_samples, test_test_samples

def make_transcript_object(samples_per_cond = 1):
    if samples_per_cond == 2:
        test_reference_samples, test_test_samples = two_sample_per_condition()
    else:
        test_reference_samples, test_test_samples = one_sample_per_condition()
    
    test_transcript_data = Transcript.Transcript(test_tx_name, test_reference_samples, test_test_samples, tx_length, max_coverage, min_coverage)
    test_transcript_data.closeAllDbs()
    return test_transcript_data

def make_txComp_results(samples_per_cond=1, weights='harmonic', context=0):
    txComp_object = TxComp.TxComp(num_reference_samples = 1,
                                  num_test_samples = 1, 
                                  methods=test_methods,
                                  sequence_context=context,
                                  sequence_context_weights=weights)

    test_transcript_data = make_transcript_object(samples_per_cond=samples_per_cond)
    total_results_dict = txComp_object.txCompare(test_transcript_data)
    return total_results_dict


def test_transcript_object():
    #test that a transcript object functions as expected
    test_reference_samples, test_test_samples = one_sample_per_condition()
    #test_reference_samples, test_test_samples = two_sample_per_condition()

    test_transcript_data = Transcript.Transcript(test_tx_name, test_reference_samples, test_test_samples, tx_length, max_coverage, min_coverage)

    assert test_transcript_data.name == 'YCR010C_mRNA'
    assert test_transcript_data.length == 852
    assert test_transcript_data._all_sample_labels == set(['WT1', 'KO1'])
    assert test_transcript_data.enoughTxCoverage(min_coverage) == True
    assert test_transcript_data.enoughTxCoverage(50) == True
    assert test_transcript_data.enoughTxCoverage(10000) == False
    assert len(test_transcript_data.getValidPositions()) == 847
    assert len(test_transcript_data.getReferenceIntensityData(pos_0)) == 89
    assert len(test_transcript_data.getReferenceIntensityData(pos_0)) == len(test_transcript_data.getReferenceDwellData(pos_0))
    assert len(test_transcript_data.getReferenceIntensityData(pos_0)) == len(test_transcript_data.getReferenceLabelData(pos_0))

    assert len(test_transcript_data.getTestIntensityData(pos_0)) == 95
    assert len(test_transcript_data.getTestIntensityData(pos_0)) == len(test_transcript_data.getTestDwellData(pos_0))
    assert len(test_transcript_data.getTestIntensityData(pos_0)) == len(test_transcript_data.getTestLabelData(pos_0))

    test_transcript_data.closeAllDbs()


def test_txComp():
    #test TxComp

    total_results_dict = pd.DataFrame(make_txComp_results()).T

    assert len(total_results_dict) == 847
    assert len(total_results_dict[total_results_dict['KS_intensity_pvalue'] <= 0.01]) == 8
    assert len(total_results_dict[total_results_dict['GMM_logit_pvalue'] <= 0.01]) == 2
    assert len(total_results_dict[abs(total_results_dict['logit_LOR']) >= 0.5]) == 94
    assert len(total_results_dict[(total_results_dict['GMM_logit_pvalue'] <= 0.01) & (abs(total_results_dict['logit_LOR']) >= 0.5)]) == 2

def exicute_SampComp_SQLDB():
    overwrite = True
    outpath = '/hps/nobackup/birney/users/logan/nanocompore_rewrite/'
    prefix = 'testing'
    results_object = results_db.SampComp_DB(outpath=outpath, prefix=prefix, overwrite=overwrite)

    total_results_dict = make_txComp_results()

    labels=set()
    for pos in total_results_dict:
        for test in total_results_dict[pos]:
            test_type = type(total_results_dict[pos][test])
            labels.add((test, test_type))

    results_object.store_test_results(test_tx_name, total_results_dict)
    results_object.closeDB()

def exicute_resultsDbManager():
    overwrite = True
    outpath = '/hps/nobackup/birney/users/logan/nanocompore_rewrite/'
    prefix = ''
    db = SampCompResultsmanager.resultsManager(outpath=outpath, prefix=prefix, overwrite=overwrite, bed_annotation=bed_fn)

    test_results_dict = test_txComp()
    db.saveData(test_tx_name, test_results_dict)

    db.finish(valid_transcripts=set(tx))

    db.closeDB()