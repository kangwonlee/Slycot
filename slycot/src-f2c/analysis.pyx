from f2c cimport *

cimport numpy as np

np.import_array()

#int sb03md_ ( char * dico, char * job, char * fact, char * trana, integer * n, doublereal * c__, integer * ldc, doublereal * a, integer * lda, doublereal * u, integer * ldu, doublereal * scale, doublereal * sep, doublereal * ferr, doublereal * wr, doublereal * wi, integer * iwork, doublereal * dwork, integer * ldwork, integer * info, ftnlen dico_len, ftnlen job_len, ftnlen fact_len, ftnlen trana_len )
cdef extern from "SB03MD.c":
    int sb03md_ ( char * dico, char * job, char * fact, char * trana, integer * n, doublereal * c__, integer * ldc, doublereal * a, integer * lda, doublereal * u, integer * ldu, doublereal * scale, doublereal * sep, doublereal * ferr, doublereal * wr, doublereal * wi, integer * iwork, doublereal * dwork, integer * ldwork, integer * info, ftnlen dico_len, ftnlen job_len, ftnlen fact_len, ftnlen trana_len )

def sb03md(dico, job, fact, trana, np.ndarray[integer, ndim=1, mode="c"] n not None, c__, ldc, a, lda, u, ldu, scale, sep, ferr, wr, wi, iwork, dwork, ldwork, info, dico_len, job_len, fact_len, trana_len):
    return sb03md_ (dico, job, fact, trana, <integer *> np.PyArray_Data(n), c__, ldc, a, lda, u, ldu, scale, sep, ferr, wr, wi, iwork, dwork, ldwork, info, dico_len, job_len, fact_len, trana_len)
