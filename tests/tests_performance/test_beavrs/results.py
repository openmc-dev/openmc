#!/usr/bin/env python

import os

# write results to file
with open('results_test.dat','w') as fh:
    fh.write(os.listdir())
