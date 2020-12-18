import pytest
import numpy as np
from loguru import logger

from nanocompore.SampCompDB import SampCompDB
from nanocompore.common import *

# set logger lever
set_logger("debug")

tol=10e-6

@pytest.mark.parametrize("pvalues", [
    (
        [0.1,0.2,0.3,0.5],
        [0.4, 0.4, 0.4, 0.5]
    ),
    (
        [0.1, 0.01, np.nan, 0.01, 0.5, 0.4, 0.01, 0.001, np.nan, np.nan, 0.01, np.nan],
        [0.13333333, 0.016, np.nan, 0.016, 0.5, 0.45714286, 0.016, 0.008, np.nan, np.nan, 0.016, np.nan]
    ),
    (
        [np.nan, np.nan, np.nan],
        [np.nan, np.nan, np.nan]
    ),
    (
        [1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1]
    )
])
def test_combine_pvalues_hou(pvalues):
    p, expected = pvalues
    corrected = list(SampCompDB._SampCompDB__multipletests_filter_nan(p))
    print(corrected)
    print(p)
    assert corrected == [pytest.approx(i, abs=tol, nan_ok=True) for i in expected]
