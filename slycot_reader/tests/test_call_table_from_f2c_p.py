import unittest

import slycot_reader.call_table_from_f2c_p as cp


class TestF2cP(unittest.TestCase):
    def setUp(self):
        self.reader = cp.F2cpReader()

    def tearDown(self):
        del self.reader

    def test_get_function_name_pattern(self):
        # input
        sample_input = 'extern int ab09ad_(char *dico, char *job, char *equil, char *ordsel, integer *n, integer *m, integer *p, integer *nr, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, doublereal *hsv, doublereal *tol, integer *iwork, doublereal *dwork, integer *ldwork, integer *iwarn, integer *info, ftnlen dico_len, ftnlen job_len, ftnlen equil_len, ftnlen ordsel_len);'

        # function under test
        p = self.reader.get_function_name_pattern()

        # check search
        results = p.search(sample_input)
        self.assertEqual('ab09ad_', results.groups()[0])

    def test_get_second_line_pattern(self):
        # input
        sample_input_output_list = [
            {'input': '''/*:ref: lsame_ 12 4 13 13 124 124 */''', 'output': 'lsame_'},
            {'input': '''/*:ref: xerbla_ 14 3 13 4 124 */''', 'output': 'xerbla_'},
            {'input': '''/*:ref: tb01id_ 14 14 13 4 4 4 7 7 4 7 4 7 4 7 4 124 */''', 'output': 'tb01id_'},
            {'input': '''/*:ref: tb01wd_ 14 16 4 4 4 7 4 7 4 7 4 7 4 7 7 7 4 4 */''', 'output': 'tb01wd_'},
            {'input': '''/*:ref: ab09ax_ 14 27 13 13 13 4 4 4 4 7 4 7 4 7 4 7 7 4 7 4 7 4 7 4 4 4 124 124 124 */''',
             'output': 'ab09ax_'},
        ]

        # function under test
        p = self.reader.get_calling_function_pattern()

        # check search
        for sample in sample_input_output_list:
            results = p.search(sample['input'])
            self.assertEqual(sample['output'], results.groups()[0])
