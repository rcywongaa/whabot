from pydrake.solvers.mathematicalprogram import MathematicalProgram, Solve
import numpy as np
from numpy import sin, cos
import math

epsilon = 1e-10
g = 9.81
l1 = 0.8
l2 = 0.32
l3 = l1
m_2 = 0.1
m_b = 4.0
m_3 = m_2
m_w = 4.0

# theta2 = np.arcsin(0.6) # ~ 36.9 degrees
theta2 = (100/180.0*np.pi)
theta1 = np.pi - theta2
theta3 = theta2
theta21 = theta2 + theta1
theta321 = theta3 + theta21
s1 = sin(theta1)
s2 = sin(theta2)
s3 = sin(theta3)
s21 = sin(theta21)
s321 = sin(theta321)
c1 = cos(theta1)
c2 = cos(theta2)
c3 = cos(theta3)
c21 = cos(theta21)
c321 = cos(theta321)

J = np.array([
    [-l1*s1 - l2*s21 - l3*s321, -l2*s21 - l3*s321, -l3*s321],
    [l1*c1 + l2*c21 + l3*c321, l2*c21 + l3*c321, l3*c321]])
print("J = " + str(J))

tau1 = -l1*c1*m_2*g + (-l1*c1 + l2/2.0)*m_b*g + (-l1*c1 + l2 + l3*c3)*m_w*g + (-l1*c1 + l2)*m_3*g

mp = MathematicalProgram()

tau23 = mp.NewContinuousVariables(2, "tau23")
tau123 = np.append(tau1, tau23)
mp.AddCost(tau23.dot(tau23))
mp.AddLinearConstraint(J.dot(tau123)[0] <= -epsilon)
mp.AddConstraint(tau23[0]**2 <= 8.0)
mp.AddConstraint(tau23[1]**2 <= 8.0)
# mp.AddLinearConstraint(J.dot(tau123)[1] >= -20)

result = Solve(mp)
print(result.is_success())
tau_result = np.append(tau1, result.GetSolution(tau23))
print("tau = " + str(tau_result))
print("F = " + str(J.dot(tau_result)))
