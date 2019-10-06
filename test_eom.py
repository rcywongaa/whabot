import unittest
import numpy as np

import eom
from constants import *

class TestEOM(unittest.TestCase):
    def test_horizontal_hold(self):
        theta1 = 0.0
        theta2 = 0.0
        theta3 = 0.0
        theta1_d = 0.0
        theta2_d = 0.0
        theta3_d = 0.0
        tau1 = l_1*m_2*g + (l_1+l_2)*m_3*g + (l_1+0.5*l_2)*m_b*g + (l_1+l_2+l_3)*m_4*g
        tau2 = l_2*m_3*g + 0.5*l_2*m_b*g + (l_2+l_3)*m_4*g
        tau3 = l_3*m_4*g
        F = 0.0
        self.assertAlmostEqual(eom.calc_theta1_dd(
            theta1, theta2, theta3,
            theta1_d, theta2_d, theta3_d,
            tau1, tau2, tau3, F,
            is_symbolic = False), 0.0)
        self.assertAlmostEqual(eom.calc_theta2_dd(
            theta1, theta2, theta3,
            theta1_d, theta2_d, theta3_d,
            tau1, tau2, tau3, F,
            is_symbolic = False), 0.0)
        self.assertAlmostEqual(eom.calc_theta3_dd(
            theta1, theta2, theta3,
            theta1_d, theta2_d, theta3_d,
            tau1, tau2, tau3, F,
            is_symbolic = False), 0.0)

if __name__ == '__main__':
    unittest.main()

