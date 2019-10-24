import unittest
import numpy as np
import pdb

import eom
from constants import *

class TestEOM(unittest.TestCase):
    def test_trapezium_configuration(self):
        theta2 = np.pi/3
        theta3 = np.pi/3
        tau2 = 0.0
        tau3 = 0.0
        F = m_4*g + (m_2 + m_3 + m_b)/2.0*g
        tau4 = F/w_r
        theta_dd = eom.solve_calc_theta_dd(
            theta2, theta3,
            0.0, 0.0,
            tau2, tau3, tau4)
        self.assertAlmostEqual(theta_dd[0], 0.0, delta=1e-3)
        self.assertAlmostEqual(theta_dd[1], 0.0, delta=1e-3)
        self.assertAlmostEqual(theta_dd[2], 0.0, delta=1e-3)

class TestFrontWheelPosition(unittest.TestCase):
    def test_zero_configuration(self):
        x, y = eom.findFrontWheelPosition(0.0, 0.0, 0.0, is_symbolic = False)
        self.assertAlmostEqual(y, 0.0)
        self.assertAlmostEqual(x, l_1 + l_2 + l_3)

    def test_theta1_90(self):
        x, y = eom.findFrontWheelPosition(np.pi/2.0, 0.0, 0.0, is_symbolic = False)
        self.assertAlmostEqual(y, -(l_1 + l_2 + l_3))
        self.assertAlmostEqual(x, 0.0)

    def test_theta2_90(self):
        x, y = eom.findFrontWheelPosition(0.0, np.pi/2.0, 0.0, is_symbolic = False)
        self.assertAlmostEqual(y, -(l_2 + l_3))
        self.assertAlmostEqual(x, l_1)

    def test_theta3_90(self):
        x, y = eom.findFrontWheelPosition(0.0, 0.0, np.pi/2.0, is_symbolic = False)
        self.assertAlmostEqual(y, -l_3)
        self.assertAlmostEqual(x, l_1 + l_2)

    def test_1(self):
        x, y = eom.findFrontWheelPosition(-3.71739739, 11.418647, -14.58469073, is_symbolic = False)
        self.assertAlmostEqual(y, -0.3)

class TestFindTheta1(unittest.TestCase):
    def test_trapezium_configuration(self):
        theta2 = np.pi/3
        theta3 = np.pi/3
        self.assertAlmostEqual(eom.findTheta1(theta2, theta3, STEP_POSITION - w_r), -np.pi/3)

if __name__ == '__main__':
    unittest.main()

