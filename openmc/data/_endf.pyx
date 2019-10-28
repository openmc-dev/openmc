# cython: c_string_type=str, c_string_encoding=ascii

cdef extern from "endf.c":
    double cfloat_endf(const char* buffer, int n)

def float_endf(s):
    cdef const char* c_string = s
    return cfloat_endf(c_string, len(s))
