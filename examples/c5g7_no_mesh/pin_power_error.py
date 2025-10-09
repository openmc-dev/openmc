import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import openmc

def reformat(arr, dev):
    arr = np.rot90(arr, 3)
    dev = np.rot90(dev, 3)
    val = arr.sum()
    factor = ((17 * 17 * 4) - (25 * 4)) / val
    arr = arr * factor
    dev = dev * factor
    return arr[0:34,0:34], dev[0:34,0:34]

def get_err(fname):
    # Load the statepoint file
    sp = openmc.StatePoint(fname)
    tally = sp.get_tally(scores=['flux'])
    flux = tally.get_slice(scores=['flux'])
    fission = tally.get_slice(scores=['fission'])

    dim = 51
    flux.std_dev.shape = (dim, dim)
    flux.mean.shape = (dim, dim)
    fission.std_dev.shape = (dim, dim)
    fission.mean.shape = (dim, dim)

    trrm_power, trrm_std_dev = reformat(fission.mean, fission.std_dev)

    np.savetxt("trrm.csv", trrm_power, delimiter=",")

    ref = np.loadtxt("reference.txt")
    ref = ref.reshape((34, 34))
    val = ref.sum()
    factor = ((17 * 17 * 4) - (25 * 4)) / val
    ref = ref * factor

    trrm_power[trrm_power == 0] = np.nan
    ref[ref == 0] = np.nan

    rel_diff = (trrm_power - ref) / ref

    rel_diff_squared = np.power(rel_diff, 2)

    RMS = math.sqrt(np.nanmean(rel_diff_squared))

    abs_percent_diff = np.absolute(rel_diff)

    print("RMS  = ", 100.0 * RMS)
    print("AAPE = ", 100.0 * np.nanmean(abs_percent_diff))
    print("MAX  = ", 100.0 *np.nanmax(abs_percent_diff))

    
n = len(sys.argv)
if (n != 2):
    print("Must provide statepoint name to read as argument.")
    exit(1)

fname = sys.argv[1]
print("Opening statepoint file: ", fname)

get_err(fname)
