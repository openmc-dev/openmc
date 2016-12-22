#!/usr/bin/env python

import os

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
        f = open(file, 'r')
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

def zaids_symbols(urr_files):
    """Return ZAIDs and isotope symbols for an ENDF/B-VII.1 ENDF-6 file list"""
    zaids = []
    symbols = []
    for f in urr_files:
        if f[len(f)-7:len(f)-5] == 'm1':
            A = f[len(f)-10:len(f)-5]
            X = f[len(f)-13:len(f)-11]
        else:
            A = f[len(f)-8:len(f)-5]
            X = f[len(f)-11:len(f)-9]
        Z = f[2:5]
        if Z[0] == '0': Z = Z[1:]
        if Z[0] == '0': Z = Z[1:]
        zero_padded_A = A
        if A[0] == '0': A = A[1:]
        if A[0] == '0': A = A[1:]
        if X[0] == '_': X = X[1:]
        symbol = X + '-' + A
        zaid = Z + zero_padded_A
        zaids.append(zaid)
        symbols.append(symbol)
    return zaids, symbols

def materials_xml(symbols):
    """Write materials.xml XML elements for URR isotopes"""
    matfile = open('materials.xml', 'w')
    for i in range(len(symbols)):
        matfile.write('<nuclide name="'+symbols[i]+'" ao="1.0" />'+'\n')
    matfile.close()
