import sys

from unittest.mock import patch

import pytest

from nanocompore.__main__ import main


def test_main_no_args_shows_help(capsys):
    with pytest.raises(SystemExit) as exception:
        with patch.object(sys, 'argv', ['nanocompore']):
            main()

    assert exception.type == SystemExit
    assert exception.value.code == 0
    captured = capsys.readouterr()
    assert "nanocompore implements the following subcommands" in captured.out

