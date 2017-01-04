#!/usr/bin/env python

import os
from subprocess import call

os.chdir('./regression')
call(['python', 'clearRegressionResults.py'])
