import numpy as np

from nanocompore.postprocessing import Postprocessor


def test_multipletests_filter_nan():
    pvals = np.array([0.1, 0.01, np.nan, 0.01, 0.5, 0.4, 0.01, 0.001, np.nan, np.nan, 0.01, np.nan])
    expected_qvals = np.array([0.13333333, 0.016, np.nan, 0.016, 0.5, 0.45714286, 0.016, 0.008, np.nan, np.nan, 0.016, np.nan])

    qvals = Postprocessor._multipletests_filter_nan(pvals)

    assert np.allclose(qvals, expected_qvals, equal_nan=True)

