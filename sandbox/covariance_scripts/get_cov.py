import os
import endf

data_dir = "/Users/fosterd/nuclear_data/ENDF-B-VIII.0_neutrons/"
files = [f for f in os.listdir(data_dir) if f.endswith(".endf") or f.endswith(".txt")]

for fname in files:
    path = os.path.join(data_dir, fname)
    mat = endf.Material(path)

    for index in mat.sections:
        mf, mt = index
        if mf == 35:  # looking for MF=35
            if (mf, mt) in mat.section_data:
                block = mat.section_data[mf, mt]
                if "subsections" in block and len(block["subsections"]) > 0:
                    print(fname, mf, mt, "has data")
                else:
                    print(fname, mf, mt, "present but empty")
            else:
                print(fname, mf, mt, "present but not parsed by endf library")
