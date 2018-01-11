from f2c cimport *

cimport numpy as np

np.import_array()

#int sb03md_ ( char * dico, char * job, char * fact, char * trana,
#       integer * n, doublereal * c__, integer * ldc, doublereal * a,
#       integer * lda, doublereal * u, integer * ldu, doublereal * scale,
#       doublereal * sep, doublereal * ferr, doublereal * wr, doublereal * wi,
#       integer * iwork, doublereal * dwork, integer * ldwork, integer * info,
#       ftnlen dico_len, ftnlen job_len, ftnlen fact_len, ftnlen trana_len )
cdef extern from "SB03MD.h":
    int sb03md_ ( char * dico, char * job, char * fact, char * trana, integer * n, doublereal * c__, integer * ldc, doublereal * a, integer * lda, doublereal * u, integer * ldu, doublereal * scale, doublereal * sep, doublereal * ferr, doublereal * wr, doublereal * wi, integer * iwork, doublereal * dwork, integer * ldwork, integer * info, ftnlen dico_len, ftnlen job_len, ftnlen fact_len, ftnlen trana_len )

def sb03md(dico, job, fact, trana,
        np.ndarray[integer, ndim=1, mode="c"] n not None,
        np.ndarray[doublereal, ndim=1, mode="c"] c__ not None,
        np.ndarray[integer, ndim=1, mode="c"] ldc not None,
        np.ndarray[doublereal, ndim=1, mode="c"] a not None,
        np.ndarray[integer, ndim=1, mode="c"] lda not None,
        np.ndarray[doublereal, ndim=1, mode="c"] u not None,
        np.ndarray[integer, ndim=1, mode="c"] ldu not None,
        np.ndarray[doublereal, ndim=1, mode="c"] scale not None,
        np.ndarray[doublereal, ndim=1, mode="c"] sep not None,
        np.ndarray[doublereal, ndim=1, mode="c"] ferr not None,
        np.ndarray[doublereal, ndim=1, mode="c"] wr not None,
        np.ndarray[doublereal, ndim=1, mode="c"] wi not None,
        np.ndarray[integer, ndim=1, mode="c"] iwork not None,
        np.ndarray[doublereal, ndim=1, mode="c"] dwork not None,
        np.ndarray[integer, ndim=1, mode="c"] ldwork not None,
        np.ndarray[integer, ndim=1, mode="c"] info,
        dico_len, job_len, fact_len, trana_len):
    sb03md_ (dico, job, fact, trana,
        <integer *> np.PyArray_Data(n),
        <doublereal *> np.PyArray_Data(c__),
        <integer *> np.PyArray_Data(ldc),
        <doublereal *> np.PyArray_Data(a),
        <integer *> np.PyArray_Data(lda),
        <doublereal *> np.PyArray_Data(u),
        <integer *> np.PyArray_Data(ldu),
        <doublereal *> np.PyArray_Data(scale),
        <doublereal *> np.PyArray_Data(sep),
        <doublereal *> np.PyArray_Data(ferr),
        <doublereal *> np.PyArray_Data(wr),
        <doublereal *> np.PyArray_Data(wi),
        <integer *> np.PyArray_Data(iwork),
        <doublereal *> np.PyArray_Data(dwork),
        <integer *> np.PyArray_Data(ldwork),
        <integer *> np.PyArray_Data(info), dico_len, job_len, fact_len, trana_len)
