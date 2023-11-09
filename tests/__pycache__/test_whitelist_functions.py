#!/usr/bin/env python 
import pytest
import nanocompore.Whitelist as Whitelist

@pytest.mark.parametrize("ref_ids, select_ref_ids, exclude_ref_ids, expected_output", [
    (['ref_id1', 'ref_id2', 'ref_id3'], [], [], ['ref_id1', 'ref_id2', 'ref_id3']),
    (['ref_id1', 'ref_id2', 'ref_id3'], ['ref_id2'], [], ['ref_id2']),
    (['ref_id1', 'ref_id2', 'ref_id3'], [], ['ref_id2'], ['ref_id1', 'ref_id3']),
    (['ref_id1', 'ref_id2', 'ref_id3'], ['ref_id1'], ['ref_id2'], ['ref_id1'])
])
def test_whitelist_private_function_cross_check(ref_ids, select_ref_ids, exclude_ref_ids, expected_output):
    result = my_function(ref_ids, select_ref_ids, exclude_ref_ids)
    assert result == expected_output