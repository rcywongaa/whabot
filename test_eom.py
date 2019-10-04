import unittest
import eom

import yaml
with open("res/constants.yaml", 'r') as stream:
    constants = yaml.safe_load(stream)

g = 9.81

m1 = constants['wheel_mass']
m2 = constants['motor1_mass']
m3 = constants['motor2_mass']
m4 = constants['wheel_mass']
mb = constants['battery_mass'] # Note we treat battery motor and battery mass as a single mass
l1 = constants['link1_length']
l2 = constants['link2_length']
l3 = constants['link3_length']
lb = constants['linkb_length']

class TestEOM(unittest.TestCase):
    def test_horizontal_hold(self):
        theta1 = 0.0
        theta2 = 0.0
        theta3 = 0.0
        theta1_d = 0.0
        theta2_d = 0.0
        theta3_d = 0.0
        tau1 = l1*m2*g + (l1+l2)*m3*g + (l1+0.5*l2)*mb*g + (l1+l2+l3)*m4*g
        tau2 = l2*m3*g + 0.5*l2*mb*g + (l2+l3)*m4*g
        tau3 = l3*m4*g
        F = 0.0
        self.assertAlmostEqual(eom.calc_theta1_dd(
            theta1, theta2, theta3,
            theta1_d, theta2_d, theta3_d,
            tau1, tau2, tau3, F), 0.0)
        self.assertAlmostEqual(eom.calc_theta2_dd(
            theta1, theta2, theta3,
            theta1_d, theta2_d, theta3_d,
            tau1, tau2, tau3, F), 0.0)
        self.assertAlmostEqual(eom.calc_theta3_dd(
            theta1, theta2, theta3,
            theta1_d, theta2_d, theta3_d,
            tau1, tau2, tau3, F), 0.0)

if __name__ == '__main__':
    unittest.main()

