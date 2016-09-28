import os
import unittest

from collections import namedtuple

from fputil import calc_a
from tests.unittests.util.fits import create_fits_with_header

class FpUtilTestCase(unittest.TestCase):
    CalcAParams = namedtuple('CalcAParams', ['ET_STATE', 'ET1B', 'ET1F', 'ET1Z', 'ET2B', 'ET2F', 'ET2Z'])

    def calc_A_fits_file_with_create(self, params):
        """Return the path to the FITS file for the given parameters.

        If the file doesn't exist, it is created.

        Params:
        -------
        params: CalcAParams
            Parameters.

        Returns:
        --------
        The path to the FITS file.

        """

        path = 'data/calc_A/{params}.fits'.format(params='_'.join(str(value).replace(' ', '').replace('&', 'and')
                                                                  for _, value in params._asdict().items()))

        if not os.path.isfile(path):
            header = {key.replace('_', '-'):value for key, value in params._asdict().items() if value is not None}
            with open(path, 'w') as f:
                create_fits_with_header(header, f)

        return path


    def test_calc_A(self):
        # use etalon 1
        self.perform_calc_A_test('S2 - Etalon 1', 7, 10000, 3, 1, 1, 1, 500, 6500, 6487.119928212755)
        self.perform_calc_A_test('S3 - Etalon 2', 1, 1, 1, 7, 10000, 3, 500, 6500, 6487.119928212755)

        # use etalon 2
        self.perform_calc_A_test('S4 - Etalon 1 & 2', 7, 10000, 3, 1, 1, 1, 500, 6500, 6487.119928212755)

        etalon_states = {'S2 - Etalon 1': 1, 'S3 - Etalon 2': 2, 'S4 - Etalon 1 & 2': 1}
        for state, etalon in etalon_states.items():
            # missing header cards for the required etalon result in an error
            params = dict(ET_STATE=state,
                          ET1B=2,
                          ET2B=4,
                          ET1F=56,
                          ET2F=43,
                          ET1Z=3,
                          ET2Z=5)
            for key in ['ET_STATE', 'ET' + str(etalon) + 'B', 'ET'+ str(etalon) + 'F', 'ET' + str(etalon) + 'Z']:
                p = {k:params[k] if k != key else None for k in params}
                pp = self.CalcAParams(**p)
                with self.assertRaises(ValueError):
                    fits_path = self.calc_A_fits_file_with_create(pp)
                    calc_a.A(fits_path, 450, 5000)

            # but there is no problem with missing cards if the respective etalon is not needed for calculating A
            other_etalon = 2 if etalon == 1 else 1
            for key in ['ET' + str(other_etalon) + 'B', 'ET' + str(other_etalon) + 'F', 'ET' + str(other_etalon) + 'Z']:
                p = {k:params[k] if k != key else None for k in params}
                pp = self.CalcAParams(**p)
                fits_path = self.calc_A_fits_file_with_create(pp)
                calc_a.A(fits_path, 450, 5000)


    def perform_calc_A_test(self, ET_STATE, ET1B, ET1F, ET1Z, ET2B, ET2F, ET2Z, R, W_line, expected_A):
        params = self.CalcAParams(ET_STATE=ET_STATE,
                                  ET1B=ET1B,
                                  ET1F=ET1F,
                                  ET1Z=ET1Z,
                                  ET2B=ET2B,
                                  ET2F=ET2F,
                                  ET2Z=ET2Z)
        fits_path = self.calc_A_fits_file_with_create(params)

        actual_A = calc_a.A(fits_path, R, W_line)
        self.assertAlmostEqual(expected_A, actual_A, places=6)

