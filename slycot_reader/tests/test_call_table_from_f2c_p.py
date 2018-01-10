import unittest

import slycot_reader.call_table_from_f2c_p as cp


class TestF2cP(unittest.TestCase):
    def setUp(self):
        self.reader = cp.F2cpReader()

    def test_get_function_name_pattern(self):
        # input
        sample_input = 'extern int ab09ad_(char *dico, char *job, char *equil, char *ordsel, integer *n, integer *m, integer *p, integer *nr, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, doublereal *hsv, doublereal *tol, integer *iwork, doublereal *dwork, integer *ldwork, integer *iwarn, integer *info, ftnlen dico_len, ftnlen job_len, ftnlen equil_len, ftnlen ordsel_len);'

        # function under test
        p = self.reader.get_function_name_pattern()

        # check search
        results = p.search(sample_input)
        self.assertEqual('ab09ad_', results.groups()[0])
