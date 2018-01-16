import unittest

import slycot_reader.write_cython as cw


class TestF2cP(unittest.TestCase):
    function_name_set = {'sb03md_', 'sb04md_', 'sg03ad_', 'sb04qd_', 'sb02md_', 'sb02mt_', 'sg02ad_', 'ab09md_',
                         'ab09md_', 'ab09nd_', 'sb10hd_', 'sb10hd_', 'sb10hd_', 'sb03od_', 'tb01pd_', 'td04ad_',
                         'td04ad_', 'sb02od_'}
    input_dict = {
        'sb03md_': {'name': 'sb03md_', 'return type': 'int', '# arg': 24,
                    'arg list': [('char *', 'dico'), ('char *', 'job'), ('char *', 'fact'), ('char *', 'trana'),
                                 ('integer *', 'n'), ('doublereal *', 'c__'), ('integer *', 'ldc'),
                                 ('doublereal *', 'a'), ('integer *', 'lda'), ('doublereal *', 'u'),
                                 ('integer *', 'ldu'), ('doublereal *', 'scale'), ('doublereal *', 'sep'),
                                 ('doublereal *', 'ferr'), ('doublereal *', 'wr'), ('doublereal *', 'wi'),
                                 ('integer *', 'iwork'), ('doublereal *', 'dwork'), ('integer *', 'ldwork'),
                                 ('integer *', 'info'), ('ftnlen', 'dico_len'), ('ftnlen', 'job_len'),
                                 ('ftnlen', 'fact_len'), ('ftnlen', 'trana_len')], 'lib': 'slycot',
                    'path': '..\\slycot\\src-f2c'},
        'sb04md_': {'name': 'sb04md_', 'return type': 'int', '# arg': 14,
                    'arg list': [('integer *', 'n'), ('integer *', 'm'), ('doublereal *', 'a'), ('integer *', 'lda'),
                                 ('doublereal *', 'b'), ('integer *', 'ldb'), ('doublereal *', 'c__'),
                                 ('integer *', 'ldc'), ('doublereal *', 'z__'), ('integer *', 'ldz'),
                                 ('integer *', 'iwork'), ('doublereal *', 'dwork'), ('integer *', 'ldwork'),
                                 ('integer *', 'info')], 'lib': 'slycot', 'path': '..\\slycot\\src-f2c'},
        'sg03ad_': {'name': 'sg03ad_', 'return type': 'int', '# arg': 31,
                    'arg list': [('char *', 'dico'), ('char *', 'job'), ('char *', 'fact'), ('char *', 'trans'),
                                 ('char *', 'uplo'), ('integer *', 'n'), ('doublereal *', 'a'), ('integer *', 'lda'),
                                 ('doublereal *', 'e'), ('integer *', 'lde'), ('doublereal *', 'q'),
                                 ('integer *', 'ldq'), ('doublereal *', 'z__'), ('integer *', 'ldz'),
                                 ('doublereal *', 'x'), ('integer *', 'ldx'), ('doublereal *', 'scale'),
                                 ('doublereal *', 'sep'), ('doublereal *', 'ferr'), ('doublereal *', 'alphar'),
                                 ('doublereal *', 'alphai'), ('doublereal *', 'beta'), ('integer *', 'iwork'),
                                 ('doublereal *', 'dwork'), ('integer *', 'ldwork'), ('integer *', 'info'),
                                 ('ftnlen', 'dico_len'), ('ftnlen', 'job_len'), ('ftnlen', 'fact_len'),
                                 ('ftnlen', 'trans_len'), ('ftnlen', 'uplo_len')], 'lib': 'slycot',
                    'path': '..\\slycot\\src-f2c'},
        'sb04qd_': {'name': 'sb04qd_', 'return type': 'int', '# arg': 14,
                    'arg list': [('integer *', 'n'), ('integer *', 'm'), ('doublereal *', 'a'), ('integer *', 'lda'),
                                 ('doublereal *', 'b'), ('integer *', 'ldb'), ('doublereal *', 'c__'),
                                 ('integer *', 'ldc'), ('doublereal *', 'z__'), ('integer *', 'ldz'),
                                 ('integer *', 'iwork'), ('doublereal *', 'dwork'), ('integer *', 'ldwork'),
                                 ('integer *', 'info')], 'lib': 'slycot', 'path': '..\\slycot\\src-f2c'},
        'sb02md_': {'name': 'sb02md_', 'return type': 'int', '# arg': 29,
                    'arg types': [13, 13, 13, 13, 13, 4, 7, 4, 7, 4, 7, 4, 7, 7, 7, 7, 4, 7, 4, 4, 7, 4, 12, 4, 124,
                                  124, 124, 124, 124], 'lib': 'slycot', 'path': '..\\slycot\\src-f2c',
                    'arg list': [('char *', 'dico'), ('char *', 'hinv'), ('char *', 'uplo'), ('char *', 'scal'),
                                 ('char *', 'sort'), ('integer *', 'n'), ('doublereal *', 'a'), ('integer *', 'lda'),
                                 ('doublereal *', 'g'), ('integer *', 'ldg'), ('doublereal *', 'q'),
                                 ('integer *', 'ldq'), ('doublereal *', 'rcond'), ('doublereal *', 'wr'),
                                 ('doublereal *', 'wi'), ('doublereal *', 's'), ('integer *', 'lds'),
                                 ('doublereal *', 'u'), ('integer *', 'ldu'), ('integer *', 'iwork'),
                                 ('doublereal *', 'dwork'), ('integer *', 'ldwork'), ('logical *', 'bwork'),
                                 ('integer *', 'info'), ('ftnlen', 'dico_len'), ('ftnlen', 'hinv_len'),
                                 ('ftnlen', 'uplo_len'), ('ftnlen', 'scal_len'), ('ftnlen', 'sort_len')]},
        'sb02mt_': {'name': 'sb02mt_', 'return type': 'int', '# arg': 28,
                    'arg types': [13, 13, 13, 13, 4, 4, 7, 4, 7, 4, 7, 4, 7, 4, 7, 4, 4, 4, 7, 4, 4, 7, 4, 4, 124, 124,
                                  124, 124], 'lib': 'slycot', 'path': '..\\slycot\\src-f2c',
                    'arg list': [('char *', 'jobg'), ('char *', 'jobl'), ('char *', 'fact'), ('char *', 'uplo'),
                                 ('integer *', 'n'), ('integer *', 'm'), ('doublereal *', 'a'), ('integer *', 'lda'),
                                 ('doublereal *', 'b'), ('integer *', 'ldb'), ('doublereal *', 'q'),
                                 ('integer *', 'ldq'), ('doublereal *', 'r__'), ('integer *', 'ldr'),
                                 ('doublereal *', 'l'), ('integer *', 'ldl'), ('integer *', 'ipiv'),
                                 ('integer *', 'oufact'), ('doublereal *', 'g'), ('integer *', 'ldg'),
                                 ('integer *', 'iwork'), ('doublereal *', 'dwork'), ('integer *', 'ldwork'),
                                 ('integer *', 'info'), ('ftnlen', 'jobg_len'), ('ftnlen', 'jobl_len'),
                                 ('ftnlen', 'fact_len'), ('ftnlen', 'uplo_len')]},
        'sg02ad_': {'name': 'sg02ad_', 'return type': 'int', '# arg': 50,
                    'arg list': [('char *', 'dico'), ('char *', 'jobb'), ('char *', 'fact'), ('char *', 'uplo'),
                                 ('char *', 'jobl'), ('char *', 'scal'), ('char *', 'sort'), ('char *', 'acc'),
                                 ('integer *', 'n'), ('integer *', 'm'), ('integer *', 'p'), ('doublereal *', 'a'),
                                 ('integer *', 'lda'), ('doublereal *', 'e'), ('integer *', 'lde'),
                                 ('doublereal *', 'b'), ('integer *', 'ldb'), ('doublereal *', 'q'),
                                 ('integer *', 'ldq'), ('doublereal *', 'r__'), ('integer *', 'ldr'),
                                 ('doublereal *', 'l'), ('integer *', 'ldl'), ('doublereal *', 'rcondu'),
                                 ('doublereal *', 'x'), ('integer *', 'ldx'), ('doublereal *', 'alfar'),
                                 ('doublereal *', 'alfai'), ('doublereal *', 'beta'), ('doublereal *', 's'),
                                 ('integer *', 'lds'), ('doublereal *', 't'), ('integer *', 'ldt'),
                                 ('doublereal *', 'u'), ('integer *', 'ldu'), ('doublereal *', 'tol'),
                                 ('integer *', 'iwork'), ('doublereal *', 'dwork'), ('integer *', 'ldwork'),
                                 ('logical *', 'bwork'), ('integer *', 'iwarn'), ('integer *', 'info'),
                                 ('ftnlen', 'dico_len'), ('ftnlen', 'jobb_len'), ('ftnlen', 'fact_len'),
                                 ('ftnlen', 'uplo_len'), ('ftnlen', 'jobl_len'), ('ftnlen', 'scal_len'),
                                 ('ftnlen', 'sort_len'), ('ftnlen', 'acc_len')], 'lib': 'slycot',
                    'path': '..\\slycot\\src-f2c'},
        'ab09md_': {'name': 'ab09md_', 'return type': 'int', '# arg': 27,
                    'arg list': [('char *', 'dico'), ('char *', 'job'), ('char *', 'equil'), ('char *', 'ordsel'),
                                 ('integer *', 'n'), ('integer *', 'm'), ('integer *', 'p'), ('integer *', 'nr'),
                                 ('doublereal *', 'alpha'), ('doublereal *', 'a'), ('integer *', 'lda'),
                                 ('doublereal *', 'b'), ('integer *', 'ldb'), ('doublereal *', 'c__'),
                                 ('integer *', 'ldc'), ('integer *', 'ns'), ('doublereal *', 'hsv'),
                                 ('doublereal *', 'tol'), ('integer *', 'iwork'), ('doublereal *', 'dwork'),
                                 ('integer *', 'ldwork'), ('integer *', 'iwarn'), ('integer *', 'info'),
                                 ('ftnlen', 'dico_len'), ('ftnlen', 'job_len'), ('ftnlen', 'equil_len'),
                                 ('ftnlen', 'ordsel_len')], 'lib': 'slycot', 'path': '..\\slycot\\src-f2c'},
        'ab09nd_': {'name': 'ab09nd_', 'return type': 'int', '# arg': 30,
                    'arg list': [('char *', 'dico'), ('char *', 'job'), ('char *', 'equil'), ('char *', 'ordsel'),
                                 ('integer *', 'n'), ('integer *', 'm'), ('integer *', 'p'), ('integer *', 'nr'),
                                 ('doublereal *', 'alpha'), ('doublereal *', 'a'), ('integer *', 'lda'),
                                 ('doublereal *', 'b'), ('integer *', 'ldb'), ('doublereal *', 'c__'),
                                 ('integer *', 'ldc'), ('doublereal *', 'd__'), ('integer *', 'ldd'),
                                 ('integer *', 'ns'), ('doublereal *', 'hsv'), ('doublereal *', 'tol1'),
                                 ('doublereal *', 'tol2'), ('integer *', 'iwork'), ('doublereal *', 'dwork'),
                                 ('integer *', 'ldwork'), ('integer *', 'iwarn'), ('integer *', 'info'),
                                 ('ftnlen', 'dico_len'), ('ftnlen', 'job_len'), ('ftnlen', 'equil_len'),
                                 ('ftnlen', 'ordsel_len')], 'lib': 'slycot', 'path': '..\\slycot\\src-f2c'},
        'sb10hd_': {'name': 'sb10hd_', 'return type': 'int', '# arg': 28,
                    'arg list': [('integer *', 'n'), ('integer *', 'm'), ('integer *', 'np'), ('integer *', 'ncon'),
                                 ('integer *', 'nmeas'), ('doublereal *', 'a'), ('integer *', 'lda'),
                                 ('doublereal *', 'b'), ('integer *', 'ldb'), ('doublereal *', 'c__'),
                                 ('integer *', 'ldc'), ('doublereal *', 'd__'), ('integer *', 'ldd'),
                                 ('doublereal *', 'ak'), ('integer *', 'ldak'), ('doublereal *', 'bk'),
                                 ('integer *', 'ldbk'), ('doublereal *', 'ck'), ('integer *', 'ldck'),
                                 ('doublereal *', 'dk'), ('integer *', 'lddk'), ('doublereal *', 'rcond'),
                                 ('doublereal *', 'tol'), ('integer *', 'iwork'), ('doublereal *', 'dwork'),
                                 ('integer *', 'ldwork'), ('logical *', 'bwork'), ('integer *', 'info')],
                    'lib': 'slycot', 'path': '..\\slycot\\src-f2c'},
        'sb03od_': {'name': 'sb03od_', 'return type': 'int', '# arg': 20,
                    'arg list': [('char *', 'dico'), ('char *', 'fact'), ('char *', 'trans'), ('integer *', 'n'),
                                 ('integer *', 'm'), ('doublereal *', 'a'), ('integer *', 'lda'), ('doublereal *', 'q'),
                                 ('integer *', 'ldq'), ('doublereal *', 'b'), ('integer *', 'ldb'),
                                 ('doublereal *', 'scale'), ('doublereal *', 'wr'), ('doublereal *', 'wi'),
                                 ('doublereal *', 'dwork'), ('integer *', 'ldwork'), ('integer *', 'info'),
                                 ('ftnlen', 'dico_len'), ('ftnlen', 'fact_len'), ('ftnlen', 'trans_len')],
                    'lib': 'slycot', 'path': '..\\slycot\\src-f2c',
                    'arg types': [13, 13, 13, 4, 4, 7, 4, 7, 4, 7, 4, 7, 7, 7, 7, 4, 4, 124, 124, 124]},
        'tb01pd_': {'name': 'tb01pd_', 'return type': 'int', '# arg': 19,
                    'arg types': [13, 13, 4, 4, 4, 7, 4, 7, 4, 7, 4, 4, 7, 4, 7, 4, 4, 124, 124], 'lib': 'slycot',
                    'path': '..\\slycot\\src-f2c',
                    'arg list': [('char *', 'job'), ('char *', 'equil'), ('integer *', 'n'), ('integer *', 'm'),
                                 ('integer *', 'p'), ('doublereal *', 'a'), ('integer *', 'lda'), ('doublereal *', 'b'),
                                 ('integer *', 'ldb'), ('doublereal *', 'c__'), ('integer *', 'ldc'),
                                 ('integer *', 'nr'), ('doublereal *', 'tol'), ('integer *', 'iwork'),
                                 ('doublereal *', 'dwork'), ('integer *', 'ldwork'), ('integer *', 'info'),
                                 ('ftnlen', 'job_len'), ('ftnlen', 'equil_len')]},
        'td04ad_': {'name': 'td04ad_', 'return type': 'int', '# arg': 24,
                    'arg types': [13, 4, 4, 4, 7, 4, 7, 4, 4, 4, 7, 4, 7, 4, 7, 4, 7, 4, 7, 4, 7, 4, 4, 124],
                    'lib': 'slycot', 'path': '..\\slycot\\src-f2c',
                    'arg list': [('char *', 'rowcol'), ('integer *', 'm'), ('integer *', 'p'), ('integer *', 'index'),
                                 ('doublereal *', 'dcoeff'), ('integer *', 'lddcoe'), ('doublereal *', 'ucoeff'),
                                 ('integer *', 'lduco1'), ('integer *', 'lduco2'), ('integer *', 'nr'),
                                 ('doublereal *', 'a'), ('integer *', 'lda'), ('doublereal *', 'b'),
                                 ('integer *', 'ldb'), ('doublereal *', 'c__'), ('integer *', 'ldc'),
                                 ('doublereal *', 'd__'), ('integer *', 'ldd'), ('doublereal *', 'tol'),
                                 ('integer *', 'iwork'), ('doublereal *', 'dwork'), ('integer *', 'ldwork'),
                                 ('integer *', 'info'), ('ftnlen', 'rowcol_len')]},
        'sb02od_': {'name': 'sb02od_', 'return type': 'int', '# arg': 43,
                    'arg list': [('char *', 'dico'), ('char *', 'jobb'), ('char *', 'fact'), ('char *', 'uplo'),
                                 ('char *', 'jobl'), ('char *', 'sort'), ('integer *', 'n'), ('integer *', 'm'),
                                 ('integer *', 'p'), ('doublereal *', 'a'), ('integer *', 'lda'), ('doublereal *', 'b'),
                                 ('integer *', 'ldb'), ('doublereal *', 'q'), ('integer *', 'ldq'),
                                 ('doublereal *', 'r__'), ('integer *', 'ldr'), ('doublereal *', 'l'),
                                 ('integer *', 'ldl'), ('doublereal *', 'rcond'), ('doublereal *', 'x'),
                                 ('integer *', 'ldx'), ('doublereal *', 'alfar'), ('doublereal *', 'alfai'),
                                 ('doublereal *', 'beta'), ('doublereal *', 's'), ('integer *', 'lds'),
                                 ('doublereal *', 't'), ('integer *', 'ldt'), ('doublereal *', 'u'),
                                 ('integer *', 'ldu'), ('doublereal *', 'tol'), ('integer *', 'iwork'),
                                 ('doublereal *', 'dwork'), ('integer *', 'ldwork'), ('logical *', 'bwork'),
                                 ('integer *', 'info'), ('ftnlen', 'dico_len'), ('ftnlen', 'jobb_len'),
                                 ('ftnlen', 'fact_len'), ('ftnlen', 'uplo_len'), ('ftnlen', 'jobl_len'),
                                 ('ftnlen', 'sort_len')], 'lib': 'slycot', 'path': '..\\slycot\\src-f2c',
                    'arg types': [13, 13, 13, 13, 13, 13, 4, 4, 4, 7, 4, 7, 4, 7, 4, 7, 4, 7, 4, 7, 7, 4, 7, 7, 7, 7, 4,
                                  7, 4, 7, 4, 7, 4, 7, 4, 12, 4, 124, 124, 124, 124, 124, 124]},
    }
    prototype_dict = {
        'sb03md_': 'int sb03md_ ( char * dico, char * job, char * fact, char * trana, integer * n, doublereal * c__, integer * ldc, doublereal * a, integer * lda, doublereal * u, integer * ldu, doublereal * scale, doublereal * sep, doublereal * ferr, doublereal * wr, doublereal * wi, integer * iwork, doublereal * dwork, integer * ldwork, integer * info, ftnlen dico_len, ftnlen job_len, ftnlen fact_len, ftnlen trana_len )',
        'sb04md_': 'int sb04md_ ( integer * n, integer * m, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, doublereal * z__, integer * ldz, integer * iwork, doublereal * dwork, integer * ldwork, integer * info )',
        'sg03ad_': 'int sg03ad_ ( char * dico, char * job, char * fact, char * trans, char * uplo, integer * n, doublereal * a, integer * lda, doublereal * e, integer * lde, doublereal * q, integer * ldq, doublereal * z__, integer * ldz, doublereal * x, integer * ldx, doublereal * scale, doublereal * sep, doublereal * ferr, doublereal * alphar, doublereal * alphai, doublereal * beta, integer * iwork, doublereal * dwork, integer * ldwork, integer * info, ftnlen dico_len, ftnlen job_len, ftnlen fact_len, ftnlen trans_len, ftnlen uplo_len )',
        'sb04qd_': 'int sb04qd_ ( integer * n, integer * m, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, doublereal * z__, integer * ldz, integer * iwork, doublereal * dwork, integer * ldwork, integer * info )',
        'sb02md_': 'int sb02md_ ( char * dico, char * hinv, char * uplo, char * scal, char * sort, integer * n, doublereal * a, integer * lda, doublereal * g, integer * ldg, doublereal * q, integer * ldq, doublereal * rcond, doublereal * wr, doublereal * wi, doublereal * s, integer * lds, doublereal * u, integer * ldu, integer * iwork, doublereal * dwork, integer * ldwork, logical * bwork, integer * info, ftnlen dico_len, ftnlen hinv_len, ftnlen uplo_len, ftnlen scal_len, ftnlen sort_len )',
        'sb02mt_': 'int sb02mt_ ( char * jobg, char * jobl, char * fact, char * uplo, integer * n, integer * m, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * q, integer * ldq, doublereal * r__, integer * ldr, doublereal * l, integer * ldl, integer * ipiv, integer * oufact, doublereal * g, integer * ldg, integer * iwork, doublereal * dwork, integer * ldwork, integer * info, ftnlen jobg_len, ftnlen jobl_len, ftnlen fact_len, ftnlen uplo_len )',
        'sg02ad_': 'int sg02ad_ ( char * dico, char * jobb, char * fact, char * uplo, char * jobl, char * scal, char * sort, char * acc, integer * n, integer * m, integer * p, doublereal * a, integer * lda, doublereal * e, integer * lde, doublereal * b, integer * ldb, doublereal * q, integer * ldq, doublereal * r__, integer * ldr, doublereal * l, integer * ldl, doublereal * rcondu, doublereal * x, integer * ldx, doublereal * alfar, doublereal * alfai, doublereal * beta, doublereal * s, integer * lds, doublereal * t, integer * ldt, doublereal * u, integer * ldu, doublereal * tol, integer * iwork, doublereal * dwork, integer * ldwork, logical * bwork, integer * iwarn, integer * info, ftnlen dico_len, ftnlen jobb_len, ftnlen fact_len, ftnlen uplo_len, ftnlen jobl_len, ftnlen scal_len, ftnlen sort_len, ftnlen acc_len )',
        'ab09md_': 'int ab09md_ ( char * dico, char * job, char * equil, char * ordsel, integer * n, integer * m, integer * p, integer * nr, doublereal * alpha, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, integer * ns, doublereal * hsv, doublereal * tol, integer * iwork, doublereal * dwork, integer * ldwork, integer * iwarn, integer * info, ftnlen dico_len, ftnlen job_len, ftnlen equil_len, ftnlen ordsel_len )',
        'ab09nd_': 'int ab09nd_ ( char * dico, char * job, char * equil, char * ordsel, integer * n, integer * m, integer * p, integer * nr, doublereal * alpha, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, doublereal * d__, integer * ldd, integer * ns, doublereal * hsv, doublereal * tol1, doublereal * tol2, integer * iwork, doublereal * dwork, integer * ldwork, integer * iwarn, integer * info, ftnlen dico_len, ftnlen job_len, ftnlen equil_len, ftnlen ordsel_len )',
        'sb10hd_': 'int sb10hd_ ( integer * n, integer * m, integer * np, integer * ncon, integer * nmeas, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, doublereal * d__, integer * ldd, doublereal * ak, integer * ldak, doublereal * bk, integer * ldbk, doublereal * ck, integer * ldck, doublereal * dk, integer * lddk, doublereal * rcond, doublereal * tol, integer * iwork, doublereal * dwork, integer * ldwork, logical * bwork, integer * info )',
        'sb03od_': 'int sb03od_ ( char * dico, char * fact, char * trans, integer * n, integer * m, doublereal * a, integer * lda, doublereal * q, integer * ldq, doublereal * b, integer * ldb, doublereal * scale, doublereal * wr, doublereal * wi, doublereal * dwork, integer * ldwork, integer * info, ftnlen dico_len, ftnlen fact_len, ftnlen trans_len )',
        'tb01pd_': 'int tb01pd_ ( char * job, char * equil, integer * n, integer * m, integer * p, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, integer * nr, doublereal * tol, integer * iwork, doublereal * dwork, integer * ldwork, integer * info, ftnlen job_len, ftnlen equil_len )',
        'td04ad_': 'int td04ad_ ( char * rowcol, integer * m, integer * p, integer * index, doublereal * dcoeff, integer * lddcoe, doublereal * ucoeff, integer * lduco1, integer * lduco2, integer * nr, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, doublereal * d__, integer * ldd, doublereal * tol, integer * iwork, doublereal * dwork, integer * ldwork, integer * info, ftnlen rowcol_len )',
        'sb02od_': 'int sb02od_ ( char * dico, char * jobb, char * fact, char * uplo, char * jobl, char * sort, integer * n, integer * m, integer * p, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * q, integer * ldq, doublereal * r__, integer * ldr, doublereal * l, integer * ldl, doublereal * rcond, doublereal * x, integer * ldx, doublereal * alfar, doublereal * alfai, doublereal * beta, doublereal * s, integer * lds, doublereal * t, integer * ldt, doublereal * u, integer * ldu, doublereal * tol, integer * iwork, doublereal * dwork, integer * ldwork, logical * bwork, integer * info, ftnlen dico_len, ftnlen jobb_len, ftnlen fact_len, ftnlen uplo_len, ftnlen jobl_len, ftnlen sort_len )',
    }
    cdef_c_func_dict = {
        'ab09nd_': 'cdef extern from "AB09ND.h":\n    int ab09nd_ ( char * dico, char * job, char * equil, char * ordsel, integer * n, integer * m, integer * p, integer * nr, doublereal * alpha, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, doublereal * d__, integer * ldd, integer * ns, doublereal * hsv, doublereal * tol1, doublereal * tol2, integer * iwork, doublereal * dwork, integer * ldwork, integer * iwarn, integer * info, ftnlen dico_len, ftnlen job_len, ftnlen equil_len, ftnlen ordsel_len )',
        'sb03od_': 'cdef extern from "SB03OD.h":\n    int sb03od_ ( char * dico, char * fact, char * trans, integer * n, integer * m, doublereal * a, integer * lda, doublereal * q, integer * ldq, doublereal * b, integer * ldb, doublereal * scale, doublereal * wr, doublereal * wi, doublereal * dwork, integer * ldwork, integer * info, ftnlen dico_len, ftnlen fact_len, ftnlen trans_len )',
        'sb02od_': 'cdef extern from "SB02OD.h":\n    int sb02od_ ( char * dico, char * jobb, char * fact, char * uplo, char * jobl, char * sort, integer * n, integer * m, integer * p, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * q, integer * ldq, doublereal * r__, integer * ldr, doublereal * l, integer * ldl, doublereal * rcond, doublereal * x, integer * ldx, doublereal * alfar, doublereal * alfai, doublereal * beta, doublereal * s, integer * lds, doublereal * t, integer * ldt, doublereal * u, integer * ldu, doublereal * tol, integer * iwork, doublereal * dwork, integer * ldwork, logical * bwork, integer * info, ftnlen dico_len, ftnlen jobb_len, ftnlen fact_len, ftnlen uplo_len, ftnlen jobl_len, ftnlen sort_len )',
        'sg03ad_': 'cdef extern from "SG03AD.h":\n    int sg03ad_ ( char * dico, char * job, char * fact, char * trans, char * uplo, integer * n, doublereal * a, integer * lda, doublereal * e, integer * lde, doublereal * q, integer * ldq, doublereal * z__, integer * ldz, doublereal * x, integer * ldx, doublereal * scale, doublereal * sep, doublereal * ferr, doublereal * alphar, doublereal * alphai, doublereal * beta, integer * iwork, doublereal * dwork, integer * ldwork, integer * info, ftnlen dico_len, ftnlen job_len, ftnlen fact_len, ftnlen trans_len, ftnlen uplo_len )',
        'sb10hd_': 'cdef extern from "SB10HD.h":\n    int sb10hd_ ( integer * n, integer * m, integer * np, integer * ncon, integer * nmeas, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, doublereal * d__, integer * ldd, doublereal * ak, integer * ldak, doublereal * bk, integer * ldbk, doublereal * ck, integer * ldck, doublereal * dk, integer * lddk, doublereal * rcond, doublereal * tol, integer * iwork, doublereal * dwork, integer * ldwork, logical * bwork, integer * info )',
        'td04ad_': 'cdef extern from "TD04AD.h":\n    int td04ad_ ( char * rowcol, integer * m, integer * p, integer * index, doublereal * dcoeff, integer * lddcoe, doublereal * ucoeff, integer * lduco1, integer * lduco2, integer * nr, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, doublereal * d__, integer * ldd, doublereal * tol, integer * iwork, doublereal * dwork, integer * ldwork, integer * info, ftnlen rowcol_len )',
        'sb04md_': 'cdef extern from "SB04MD.h":\n    int sb04md_ ( integer * n, integer * m, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, doublereal * z__, integer * ldz, integer * iwork, doublereal * dwork, integer * ldwork, integer * info )',
        'tb01pd_': 'cdef extern from "TB01PD.h":\n    int tb01pd_ ( char * job, char * equil, integer * n, integer * m, integer * p, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, integer * nr, doublereal * tol, integer * iwork, doublereal * dwork, integer * ldwork, integer * info, ftnlen job_len, ftnlen equil_len )',
        'ab09md_': 'cdef extern from "AB09MD.h":\n    int ab09md_ ( char * dico, char * job, char * equil, char * ordsel, integer * n, integer * m, integer * p, integer * nr, doublereal * alpha, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, integer * ns, doublereal * hsv, doublereal * tol, integer * iwork, doublereal * dwork, integer * ldwork, integer * iwarn, integer * info, ftnlen dico_len, ftnlen job_len, ftnlen equil_len, ftnlen ordsel_len )',
        'sg02ad_': 'cdef extern from "SG02AD.h":\n    int sg02ad_ ( char * dico, char * jobb, char * fact, char * uplo, char * jobl, char * scal, char * sort, char * acc, integer * n, integer * m, integer * p, doublereal * a, integer * lda, doublereal * e, integer * lde, doublereal * b, integer * ldb, doublereal * q, integer * ldq, doublereal * r__, integer * ldr, doublereal * l, integer * ldl, doublereal * rcondu, doublereal * x, integer * ldx, doublereal * alfar, doublereal * alfai, doublereal * beta, doublereal * s, integer * lds, doublereal * t, integer * ldt, doublereal * u, integer * ldu, doublereal * tol, integer * iwork, doublereal * dwork, integer * ldwork, logical * bwork, integer * iwarn, integer * info, ftnlen dico_len, ftnlen jobb_len, ftnlen fact_len, ftnlen uplo_len, ftnlen jobl_len, ftnlen scal_len, ftnlen sort_len, ftnlen acc_len )',
        'sb04qd_': 'cdef extern from "SB04QD.h":\n    int sb04qd_ ( integer * n, integer * m, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * c__, integer * ldc, doublereal * z__, integer * ldz, integer * iwork, doublereal * dwork, integer * ldwork, integer * info )',
        'sb03md_': 'cdef extern from "SB03MD.h":\n    int sb03md_ ( char * dico, char * job, char * fact, char * trana, integer * n, doublereal * c__, integer * ldc, doublereal * a, integer * lda, doublereal * u, integer * ldu, doublereal * scale, doublereal * sep, doublereal * ferr, doublereal * wr, doublereal * wi, integer * iwork, doublereal * dwork, integer * ldwork, integer * info, ftnlen dico_len, ftnlen job_len, ftnlen fact_len, ftnlen trana_len )',
        'sb02md_': 'cdef extern from "SB02MD.h":\n    int sb02md_ ( char * dico, char * hinv, char * uplo, char * scal, char * sort, integer * n, doublereal * a, integer * lda, doublereal * g, integer * ldg, doublereal * q, integer * ldq, doublereal * rcond, doublereal * wr, doublereal * wi, doublereal * s, integer * lds, doublereal * u, integer * ldu, integer * iwork, doublereal * dwork, integer * ldwork, logical * bwork, integer * info, ftnlen dico_len, ftnlen hinv_len, ftnlen uplo_len, ftnlen scal_len, ftnlen sort_len )',
        'sb02mt_': 'cdef extern from "SB02MT.h":\n    int sb02mt_ ( char * jobg, char * jobl, char * fact, char * uplo, integer * n, integer * m, doublereal * a, integer * lda, doublereal * b, integer * ldb, doublereal * q, integer * ldq, doublereal * r__, integer * ldr, doublereal * l, integer * ldl, integer * ipiv, integer * oufact, doublereal * g, integer * ldg, integer * iwork, doublereal * dwork, integer * ldwork, integer * info, ftnlen jobg_len, ftnlen jobl_len, ftnlen fact_len, ftnlen uplo_len )',
    }
    c_file_name_dict = {
        'sb04md_': 'SB04MD',
        'ab09md_': 'AB09MD',
        'sb03md_': 'SB03MD',
        'tb01pd_': 'TB01PD',
        'sb03od_': 'SB03OD',
        'td04ad_': 'TD04AD',
        'sb10hd_': 'SB10HD',
        'ab09nd_': 'AB09ND',
        'sb04qd_': 'SB04QD',
        'sb02od_': 'SB02OD',
        'sg02ad_': 'SG02AD',
        'sb02mt_': 'SB02MT',
        'sb02md_': 'SB02MD',
        'sg03ad_': 'SG03AD',
    }
    py_func_name_dict = {
        'sb02od_': 'sb02od',
        'ab09md_': 'ab09md',
        'sb03md_': 'sb03md',
        'sb02mt_': 'sb02mt',
        'sb03od_': 'sb03od',
        'sg03ad_': 'sg03ad',
        'sb02md_': 'sb02md',
        'sg02ad_': 'sg02ad',
        'ab09nd_': 'ab09nd',
        'tb01pd_': 'tb01pd',
        'td04ad_': 'td04ad',
        'sb10hd_': 'sb10hd',
        'sb04md_': 'sb04md',
        'sb04qd_': 'sb04qd',
    }

    def setUp(self):
        self.writer = cw.Dict2Cython(self.input_dict, self.function_name_set)

    def test_get_function_prototype_text(self):
        for function_name in self.function_name_set:
            self.assertEqual(
                self.prototype_dict[function_name],
                self.writer.get_function_prototype_text(function_name)
            )

    def test_get_cdef_c_func_block(self):
        for function_name in self.function_name_set:
            self.assertEqual(
                self.cdef_c_func_dict[function_name],
                self.writer.get_cdef_c_func_block(function_name)
            )

    def test_get_c_file_name(self):
        for function_name in self.function_name_set:
            # print('%r:%r,' % (function_name, self.writer.get_c_file_name(function_name)))
            self.assertEqual(
                self.c_file_name_dict[function_name],
                self.writer.get_c_file_name(function_name)
            )

    def test_get_py_func_name(self):
        for function_name in self.function_name_set:
            # print('%r:%r,' % (function_name, self.writer.get_py_func_name(function_name)))
            self.assertEqual(
                self.py_func_name_dict[function_name],
                self.writer.get_py_func_name(function_name)
            )


if __name__ == '__main__':
    unittest.main()
