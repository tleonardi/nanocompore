import pytest
from nanocompore.TxComp import *
from scipy.stats import combine_pvalues
import numpy as np
from unittest import mock


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
    ( (1,4,3,2,9), (1,7,8,8,5), (0.4619, 0.3277, 0.3291))
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
                                          'rep1':{
                                             'intensity': np.random.normal(loc=100, scale=10.0, size=100),
                                             'dwell': np.random.normal(loc=100, scale=10.0, size=100),
                                             'coverage': 100
                                          },
                                          'rep2':{
                                             'intensity':np.random.normal(loc=100, scale=10.0, size=100),
                                             'dwell':np.random.normal(loc=100, scale=10.0, size=100),
                                             'coverage': 100
                                          }
                                      },
                                      'KD':{
                                          'rep1':{
                                             'intensity':np.random.normal(loc=120, scale=10.0, size=100),
                                             'dwell':np.random.normal(loc=120, scale=10.0, size=100),
                                             'coverage': 100
                                          },
                                          'rep2':{
                                             'intensity':np.random.normal(loc=120, scale=10.0, size=100),
                                             'dwell':np.random.normal(loc=120, scale=10.0, size=100),
                                             'coverage': 100
                                          }
                                      }
                                  }
                                 }
    expected = {'GMM_anova': [0.009294583692737432, 0.0007186082196801213, 0.0073120479947770415, 0.0006747727949872362, 0.00607648353543618, 0.004730678594686069, 0.004780377368257902, 0.0032067424187288665, 0.00040586223756605407, 0.0009465280405188231],
		'GMM_logit': [1.2742453287653416e-39, 3.3968653938213694e-40, 1.9321679678622595e-36, 8.482777798354296e-40, 6.928304506982671e-39, 7.065038671811672e-40, 1.8392720921153275e-40, 4.6826664356268694e-32, 5.922884891638699e-34, 3.1972432623454785e-40]
		}
    return((test_ref_pos_list, expected))


def test_txComp_GMM_anova(test_ref_pos_list):
    ml = mock.Mock()
    tol=0.000000001
    res = txCompare(test_ref_pos_list[0], methods=['GMM'], logit=False, sequence_context=2, min_coverage=3, logger=ml, strict=True)
    GMM_pvalues = [pos['txComp']['GMM_pvalue'] for pos in res ]
    assert GMM_pvalues == [pytest.approx(i, abs=tol) for i in test_ref_pos_list[1]['GMM_anova']] 

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
    ml = mock.Mock()
    tol=0.000000001
    with pytest.raises(NanocomporeError):
        txCompare(test_ref_pos_list_0_var, methods=['GMM'], logit=False, sequence_context=2, min_coverage=3, logger=ml, strict=True)

def test_txComp_GMM_logit(test_ref_pos_list):
    ml = mock.Mock()
    tol=0.000000001
    res = txCompare(test_ref_pos_list[0], methods=['GMM'], logit=True, anova=False, sequence_context=2, min_coverage=3, logger=ml, strict=True)
    GMM_logit = [pos['txComp']['GMM_logit_pvalue'] for pos in res ]
    assert GMM_logit == [pytest.approx(i, abs=tol) for i in test_ref_pos_list[1]['GMM_logit']] 
