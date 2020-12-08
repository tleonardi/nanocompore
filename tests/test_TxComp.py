import sys

import numpy as np
import pytest
from scipy.stats import combine_pvalues
from loguru import logger

from nanocompore.TxComp import *
from nanocompore.common import *

# set logger lever
set_logger("debug")


@pytest.mark.parametrize("pvalues", [
    np.array([0.1,0.2,0.3,0.5]),
    np.array([0.1,0.7,0.001,0.5]),
    np.array([0.1,0.5,np.nan,0.5]),
    np.array([1,1,0.001,0.5]),
])
def test_combine_pvalues_hou(pvalues):
    weights = np.array([1,1,1,1])
    cor_mat = np.zeros((4,4))
    hou = combine_pvalues_hou(pvalues, weights, cor_mat)
    fisher = combine_pvalues(pvalues, method='fisher')[1]
    if np.isnan(hou):
        assert np.isnan(fisher)
    else:
        assert hou == fisher

@pytest.mark.parametrize("pvalues", [
    np.array([0.1,np.inf,0.3,0.5]),
    np.array([1,1,0.001,0]),
    np.array([1,7,0.001,0.1]),
])
def test_combine_pvalues_raises_exception_with_invalid_pvalues(pvalues):
    weights = np.array([1,1,1,1])
    cor_mat = np.zeros((4,4))
    with pytest.raises(NanocomporeError):
        combine_pvalues_hou(pvalues, weights, cor_mat)


@pytest.mark.parametrize("v1, v2, expected", [
    ( (1,4,3,2,9), (1,7,8,8,5), (0.4619, 0.3277, 0.3571))
])
def test_nonparametric_test(v1, v2, expected):
    """ Test the non parametric methods
        expected: a tuple with the expected results for MW, KS and TT
    """
    tol=0.0001
    assert nonparametric_test(v1, v2, v1, v2, method="MW") == (pytest.approx(expected[0], abs=tol), pytest.approx(expected[0], abs=tol))
    assert nonparametric_test(v1, v2, v1, v2, method="TT") == (pytest.approx(expected[1], abs=tol), pytest.approx(expected[1], abs=tol))
    assert nonparametric_test(v1, v2, v1, v2, method="KS") == (pytest.approx(expected[2], abs=tol), pytest.approx(expected[2], abs=tol))

@pytest.mark.parametrize("x, expected", [
        ( (1, 2, 7, 33, 10), 1243 ),
        ((2), 4 ),
        ( (10, -10, 0, 1), 201)
])
def test_sum_of_squares(x, expected):
    assert sum_of_squares(x) == expected

@pytest.fixture
def test_ref_pos_list():
    test_ref_pos_list=[None]*10
    np.random.seed(seed=6354565)
    for pos in range(0,10):
        test_ref_pos_list[pos] = {'data':{
                                      'WT':{
                                          'WT1':{
                                             'intensity': np.random.normal(loc=100, scale=10.0, size=100),
                                             'dwell': np.random.normal(loc=100, scale=10.0, size=100),
                                             'coverage': 100
                                          },
                                          'WT2':{
                                             'intensity':np.random.normal(loc=100, scale=10.0, size=100),
                                             'dwell':np.random.normal(loc=100, scale=10.0, size=100),
                                             'coverage': 100
                                          }
                                      },
                                      'KD':{
                                          'KD1':{
                                             'intensity':np.random.normal(loc=120, scale=10.0, size=100),
                                             'dwell':np.random.normal(loc=120, scale=10.0, size=100),
                                             'coverage': 100
                                          },
                                          'KD2':{
                                             'intensity':np.random.normal(loc=120, scale=10.0, size=100),
                                             'dwell':np.random.normal(loc=120, scale=10.0, size=100),
                                             'coverage': 100
                                          }
                                      }
                                  }
                                 }
    # These expected values have been checked against the R implementations of aov() and glm(family="binomial")
    expected = {'GMM_anova': [0.0008574768473501677, 0.0036329291397528157, 0.007312047981252302, 0.0017335646468135102, np.nan, 0.004197519576768562, 0.004730678586860965, 0.0028228474915020945, 0.0023262178697710987, 0.00020764199465021126],
		'GMM_logit': [1.274245328765287e-39, 3.3968653938213694e-40, 1.9321679678623975e-36, 8.482777798353687e-40, np.nan, 7.06503867181238e-40, 1.839272092115274e-40, 9.162002495725215e-32, 5.922884891638699e-34, 3.1972432623454785e-40]
		}
    return((test_ref_pos_list, expected))


def test_txComp_GMM_anova(test_ref_pos_list):
    if sys.version_info < (3, 6):
        tol = 0.001
    else:
        tol=0.00000001
    res = txCompare("ref_id", test_ref_pos_list[0], methods=['GMM'], logit=False, sequence_context=2, min_coverage=3, allow_warnings=False, random_state=np.random.RandomState(seed=42))
    GMM_pvalues = [pos['txComp']['GMM_anova_pvalue'] for pos in res ]
    assert GMM_pvalues == [pytest.approx(i, abs=tol, nan_ok=True) for i in test_ref_pos_list[1]['GMM_anova']]

def test_txComp_GMM_logit(test_ref_pos_list):
    tol=0.000000001
    res = txCompare("ref_id", test_ref_pos_list[0], methods=['GMM'], logit=True, anova=False, sequence_context=2, min_coverage=3, allow_warnings=False, random_state=np.random.RandomState(seed=42))
    GMM_logit = [pos['txComp']['GMM_logit_pvalue'] for pos in res ]
    print(GMM_logit)
    print(test_ref_pos_list[1]['GMM_logit'])
    print([pos['txComp']['GMM_model']['cluster_counts'] for pos in res ])

    assert GMM_logit == [pytest.approx(i, abs=tol, nan_ok=True) for i in test_ref_pos_list[1]['GMM_logit']]

@pytest.fixture
def test_ref_pos_list_0_var():
    test_ref_pos_list=[None]*10
    np.random.seed(seed=6354565)
    for pos in range(0,10):
        int_1 = np.random.normal(loc=100, scale=10.0, size=100)
        dwell_1 = np.random.normal(loc=100, scale=10.0, size=100)
        int_2 = np.random.normal(loc=500, scale=10.0, size=100)
        dwell_2 = np.random.normal(loc=500, scale=10.0, size=100)
        test_ref_pos_list[pos] = {'data':{
                                      'WT':{
                                          'WT1':{
                                             'intensity': int_1,
                                             'dwell': dwell_1,
                                             'coverage': 100
                                          },
                                          'WT2':{
                                             'intensity': int_1,
                                             'dwell': dwell_1,
                                             'coverage': 100
                                          }
                                      },
                                      'KD':{
                                          'KD1':{
                                             'intensity': int_2,
                                             'dwell': dwell_2,
                                             'coverage': 100
                                          },
                                          'KD2':{
                                             'intensity': int_2,
                                             'dwell': dwell_2,
                                             'coverage': 100
                                          }
                                      }
                                  }
                                 }
    return(test_ref_pos_list)

def test_txComp_GMM_anova_0_var(test_ref_pos_list_0_var):
    with pytest.raises(NanocomporeError):
        txCompare("ref_id", test_ref_pos_list_0_var, methods=['GMM'], logit=False, sequence_context=2, min_coverage=3, allow_warnings=False, random_state=np.random.RandomState(seed=42))

@pytest.fixture
def test_ref_pos_list_dup_lab():
    test_ref_pos_list = [None]
    test_ref_pos_list[0] = {'data':{
                                  'WT':{
                                      'Rep1':{
                                         'intensity': [0,0],
                                         'dwell': [0,0],
                                         'coverage': 100
                                      },
                                      'Rep2':{
                                         'intensity': [0,0],
                                         'dwell': [0,0],
                                         'coverage': 100
                                      }
                                  },
                                  'KD':{
                                      'Rep1':{
                                         'intensity': [0,0],
                                         'dwell': [0,0],
                                         'coverage': 100
                                      },
                                      'Rep2':{
                                         'intensity': [0,0],
                                         'dwell': [0,0],
                                         'coverage': 100
                                      }
                                  }
                              }
                             }
    return(test_ref_pos_list)

def test_txComp_GMM_dup_lab(test_ref_pos_list_dup_lab):
    with pytest.raises(NanocomporeError):
        txCompare("ref_id", test_ref_pos_list_dup_lab, methods=['GMM'], logit=False, sequence_context=2, min_coverage=3, allow_warnings=False, random_state=np.random.RandomState(seed=42))

def test_txComp_lowCov(test_ref_pos_list):
    """ This test ensures that txCompare runs also when the number of covered positions
        in a reference is below the threshold
    """
    test_ref_pos_list = test_ref_pos_list[0]
    low_cov_positions = [0,1,5]
    for pos in low_cov_positions:
        test_ref_pos_list[pos]['data']['WT']['WT1']['coverage'] = 1
    results = txCompare("ref_id", test_ref_pos_list, methods=['GMM'], logit=False, sequence_context=2, min_coverage=30, allow_warnings=False, random_state=np.random.RandomState(seed=42))
    for pos in results:
        if 'txComp' in pos:
            # If the original p-value was nan, the context p-value also has to be nan
            if np.isnan(pos['txComp']['GMM_anova_pvalue']):
                assert np.isnan(pos['txComp']['GMM_anova_pvalue_context_2'])
            else:
                assert not np.isnan(pos['txComp']['GMM_anova_pvalue_context_2'])
