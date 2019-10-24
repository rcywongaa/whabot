import unittest
import numpy as np
import pdb

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
        tau1 = -l_1*m_2*g - (l_1+l_2)*m_3*g - (l_1+0.5*l_2)*m_b*g - (l_1+l_2+l_3)*m_4*g
        tau2 = -l_2*m_3*g - 0.5*l_2*m_b*g - (l_2+l_3)*m_4*g
        tau3 = -l_3*m_4*g
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

    def test_trapezium_configuration(self):
        # Offset actual angle slightly...
        # because I don't want to solve the simultaneous equation
        # for required_angle and required_penetration
        required_F_x = (m_2 + m_3 + m_b)/2.0*g/np.tan(np.pi/3 - 0.0002)
        required_penetration = required_F_x / CONTACT_SPRING_STIFFNESS
        required_angle = np.arccos((required_penetration + STEP_POSITION - w_r - l_2)/(l_1+l_3))
        theta = [-required_angle, required_angle, required_angle]
        theta_d = [0., 0., 0., 0.]
        tau1 = 0.0
        tau2 = 0.0
        tau3 = 0.0
        F = m_4*g + (m_2 + m_3 + m_b)/2.0*g
        self.assertAlmostEqual(eom.calc_theta1_dd(
            theta[0], theta[1], theta[2],
            theta_d[0], theta_d[1], theta_d[2],
            tau1, tau2, tau3, F,
            is_symbolic = False), 0.0, delta=1e-3)
        self.assertAlmostEqual(eom.calc_theta2_dd(
            theta[0], theta[1], theta[2],
            theta_d[0], theta_d[1], theta_d[2],
            tau1, tau2, tau3, F,
            is_symbolic = False), 0.0, delta=1e-3)
        self.assertAlmostEqual(eom.calc_theta3_dd(
            theta[0], theta[1], theta[2],
            theta_d[0], theta_d[1], theta_d[2],
            tau1, tau2, tau3, F,
            is_symbolic = False), 0.0, delta=1e-3)

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
        theta4 = 0.0
        pdb.set_trace()
        self.assertAlmostEqual(eom.findTheta1(theta2, theta3, theta4), -np.pi/3)

if __name__ == '__main__':
    unittest.main()

