import os.path
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlopen

_BLOCK_SIZE = 16384


def download(url):
    """Download file from a URL

    Parameters
    ----------
    url : str
        URL from which to download

    Returns
    -------
    basename : str
        Name of file written locally

    """
    req = urlopen(url)

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
    return basename
