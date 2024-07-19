#!/usr/bin/env python3

# helper script that launches the openmc binary

import os
import sys
import sysconfig
from pathlib import Path


def main():
    os.execv(
        Path(sysconfig.get_path("platlib")) / "openmc" / "bin" / "openmc", sys.argv
    )


if __name__ == "__main__":
    main()
