from unittest import TestCase, skip
import warnings
from boxmodel_functions import condensation_solver_euler as condensation_solver
from boxmodel_functions import condensation
from boxmodel_functions import radius

class _TestWarnings(TestCase):
    def setUp(self):
        warnings.filterwarnings('error')

    def tearDown(self):
        warnings.resetwarnings()

class TestCondensationSolver(_TestWarnings):
    def test_negative_saturation(self):
        condensation_solver(r_old=1.02113775484e-06,
                                    dt=2,
                                    E=15.63245753,
                                    es=1599.80785594,
                                    T=287.198552482,
                                    S=-0.0115141998358)

    def test_negative_saturation_with_small_timestep(self):
        condensation_solver(r_old=1.02113775484e-06,
                                    dt=0.1,
                                    E=15.63245753,
                                    es=1599.80785594,
                                    T=287.198552482,
                                    S=-0.0115141998358)

    def test_negative_saturation_and_nonzero_qc(self):
        condensation_solver(r_old=1.02113775484e-06,
                                    dt=2,
                                    E=15.63245753,
                                    es=1599.80785594,
                                    T=287.198552482,
                                    S=-0.0115141998358)

    def test_positive_saturation(self):
        condensation_solver(r_old=1.02113775484e-06,
                                    dt=2,
                                    E=15.63245753,
                                    es=1599.80785594,
                                    T=287.198552482,
                                    S=0.0115141998358)

    def test_negative_saturation_zero_rhs(self):
        E = 16.
        r = 1.e-6
        S = -E * r
        condensation_solver(r_old=r,
                            dt=2,
                            E=E,
                            es=1599.80785594,
                            T=287.198552482,
                            S=S)

    def test_negative_saturation_positive_rhs(self):
        E = 16.
        r = 1.e-6
        S = -E * r / 2
        condensation_solver(r_old=r,
                            dt=2,
                            E=E,
                            es=1599.80785594,
                            T=287.198552482,
                            S=S)
        
    def test_negative_saturation_negative_rhs(self):
        E = 16.
        r = 1.e-6
        S = -E * r * 2
        condensation_solver(r_old=r,
                            dt=2,
                            E=E,
                            es=1599.80785594,
                            T=287.198552482,
                            S=S)

class TestRadiusCalculator(TestCase):
    def test_qc_greater_then_zero(self):
        inp = {'qc': 10.e-3 , 'N': 0.1e6}
        out = 0
        self.assertGreaterEqual(radius(**inp), out)

    def test_qc_equal_to_zero(self):
        inp = {'qc': 0. , 'N': 0.1e6}
        out = 0.
        self.assertEqual(radius(**inp), out)

    def test_qc_lower_then_zero(self):
        with self.assertRaises(ValueError):
            inp = {'qc': -10.e-3 , 'N': 0.1e6}
            radius(**inp)

    def test_qc_equal_to_none(self):
        with self.assertRaises(TypeError):
            inp = {'qc': None, 'N': 0.1e6}
            radius(**inp)

    def test_N_equal_to_zero(self):
        with self.assertRaises(ZeroDivisionError):
            inp = {'qc': 10.e-3 , 'N': 0.}
            radius(**inp)

    def test_N_equal_to_none(self):
        with self.assertRaises(TypeError):
            inp = {'qc': 10.e-3, 'N': None}
            radius(**inp)
