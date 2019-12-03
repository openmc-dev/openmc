import hashlib
import os.path
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlopen, Request

_BLOCK_SIZE = 16384


def download(url, checksum=None, as_browser=False, **kwargs):
    """Download file from a URL

    Parameters
    ----------
    url : str
        URL from which to download
    checksum : str or None
        MD5 checksum to check against
    as_browser : bool
        Change User-Agent header to appear as a browser
    kwargs : dict
        Keyword arguments passed to :func:urllib.request.urlopen

    Returns
    -------
    basename : str
        Name of file written locally

    """
    if as_browser:
        page = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    else:
        page = url
    req = urlopen(page, **kwargs)
    # Get file size from header
    file_size = req.length

    # Check if file already downloaded
    basename = Path(urlparse(url).path).name
    if os.path.exists(basename):
        if os.path.getsize(basename) == file_size:
            print('Skipping {}, already downloaded'.format(basename))
            return basename

    # Copy file to disk in chunks
    print('Downloading {}... '.format(basename), end='')
    downloaded = 0
    with open(basename, 'wb') as fh:
        while True:
            chunk = req.read(_BLOCK_SIZE)
            if not chunk:
                break
            fh.write(chunk)
            downloaded += len(chunk)
            status = '{:10}  [{:3.2f}%]'.format(
                downloaded, downloaded * 100. / file_size)
            print(status + '\b'*len(status), end='')
        print('')

    if checksum is not None:
        downloadsum = hashlib.md5(open(basename, 'rb').read()).hexdigest()
        if downloadsum != checksum:
            raise IOError("MD5 checksum for {} does not match. If this is your first "
                          "time receiving this message, please re-run the script. "
                          "Otherwise, please contact OpenMC developers by emailing "
                          "openmc-users@googlegroups.com.".format(basename))

    return basename
