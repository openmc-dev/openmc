import os
import filecmp

from openmc import _utils
import pytest

@pytest.fixture()
def download_photos(run_in_tmpdir):
    """use _utils download() function to download the same picture three times,
       twice to get unique names, & a third time to use the already downloaded
       block of code"""
    _utils.download("https://i.ibb.co/HhKFc8x/small.jpg")
    _utils.download("https://tinyurl.com/y4t38ugb")
    _utils.download("https://tinyurl.com/y4t38ugb", as_browser=True)


def test_checksum_error(run_in_tmpdir):
    """use download() in such a way that will test the checksum error line"""
    phrase = "MD5 checksum for y4t38ugb"
    with pytest.raises(OSError, match=phrase):
        _utils.download("https://tinyurl.com/y4t38ugb", as_browser=True,
                        checksum="not none")


def test_photos(download_photos):
    assert filecmp.cmp("small.jpg", "y4t38ugb")
