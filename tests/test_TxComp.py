import pytest
from nanocompore.TxComp import *
from scipy.stats import combine_pvalues
import numpy as np

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
