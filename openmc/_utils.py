import hashlib
import os.path
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlopen

_BLOCK_SIZE = 16384


def download(url, checksum=None, context=None):
    """Download file from a URL

    Parameters
    ----------
    url : str
        URL from which to download
    checksum : str or None
        MD5 checksum to check against
    context : ssl.SSLContext instance or None
        For example ssl._create_unverified_context()

    Returns
    -------
    basename : str
        Name of file written locally

    """
    req = urlopen(url, context=context)

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
