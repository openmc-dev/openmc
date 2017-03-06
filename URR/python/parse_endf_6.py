#!/usr/bin/env python

import os
import njoy as njoy

def urr_filenames(endf_6_path):
    """
    Return names of ENDF-6 files having a URR evaluation

    arguments:
    endf_6_path -- path to ENDF-6 files
    """
    urr_files = []
    EL = []
    EH = []
    for file in os.listdir(endf_6_path):
        f = open(endf_6_path+file, 'r')
        cnt = 0
        for l in f:
            fields = l.split()
            if (len(fields) > 1 and fields[len(fields)-2] == '2151'):
                cnt += 1
                if (cnt >= 3
                    and fields[2] == '2'
                    and (fields[3] == '1' or fields[3] =='2')
                    and int(fields[4]) != 6*int(fields[5])):
                    if not file in urr_files:
                        EL.append(fields[0])
                        EH.append(fields[1])
                        urr_files.append(file)
                    break
        f.close()
    return urr_files

def zaids_symbols(urr_files, evaluation):
    """Return ZAIDs and isotope symbols for an ENDF/B-VII.1 ENDF-6 file list"""
    zaids = []
    symbols = []
    for f in urr_files:
        if evaluation == 'ENDF':
            if f[len(f)-7:len(f)-5] == 'm1':
                zero_padded_A = f[len(f)-10:len(f)-6]
                X = f[len(f)-13:len(f)-11]
            else:
                zero_padded_A = f[len(f)-8:len(f)-5]
                X = f[len(f)-11:len(f)-9]
            Z = f[2:5]
            if Z[0] == '0': Z = Z[1:]
            if Z[0] == '0': Z = Z[1:]
            if X[0] == '_': X = X[1:]
        elif evaluation == 'JENDL':
            f = f[:-4]
            if f[-1] == 'm':
                zero_padded_A = f[-4:]
                X = f[0:-4]
            else:
                zero_padded_A = f[-3:]
                X = f[0:-3]
            Z = njoy.X_Z[X]
        else:
            raise ValueError('Unrecognized evaluated nuclear data library')
        A = zero_padded_A
        if A[0] == '0': A = A[1:]
        if A[0] == '0': A = A[1:]
        zaid = Z + zero_padded_A
        zaids.append(zaid)
        symbol = X + '-' + A
        symbols.append(symbol)
    return zaids, symbols
