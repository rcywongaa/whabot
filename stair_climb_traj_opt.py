from pydrake.solvers.mathematicalprogram import MathematicalProgram, Solve
import numpy as np
from numpy import sin, cos
import math

HUB_MOTOR_MAX_TORQUE = 20
MIN_NORMAL_REACTION = 4
epsilon = 1e-10
g = 9.81
l1 = 0.8
l2 = 0.32
l3 = l1
m_2 = 0.1
m_b = 4.0
m_3 = m_2
m_w = 4.0

def findMaxForce(theta1, theta2, theta3):
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

    tau1 = -l1*c1*m_2*g + (-l1*c1 + l2/2.0)*m_b*g + (-l1*c1 + l2 + l3*c3)*m_w*g + (-l1*c1 + l2)*m_3*g

    mp = MathematicalProgram()

    tau23 = mp.NewContinuousVariables(2, "tau23")
    tau123 = np.append(tau1, tau23)
    mp.AddCost(J.dot(tau123)[0])
    mp.AddConstraint(tau23[0]**2 <= 8.0)
    mp.AddConstraint(tau23[1]**2 <= 8.0)
    mp.AddLinearConstraint(J.dot(tau123)[1] >= epsilon)

    result = Solve(mp)
    is_success = result.is_success()
    torques = result.GetSolution(tau23)
    full_tau = np.append(tau1, torques)
    output_force = J.dot(full_tau)
    # print("is_success = " + str(is_success))
    # print("tau = " + str(full_tau))
    # print("F = " + str(J.dot(full_tau)))
    return is_success, output_force, torques

def solve(y):
    valid_thetas = []
    # theta2 = np.arcsin(0.6) # ~ 36.9 degrees
    for theta2 in np.arange(0.0, np.pi/2.0, 0.02): # roughly 1 degree intervals
        for theta3 in np.arange(0.0, np.pi/2.0, 0.02):
            s2 = sin(theta2)
            s3 = sin(theta3)
            s23 = sin(theta2+theta3)
            c2 = cos(theta2)
            c3 = cos(theta3)
            c23 = cos(theta2+theta3)
            # Given theta2, theta3 and y, find theta1
            # i.e. solve                                          l1*s1 + l2*s12 + l3*s123 = y
            # Use sum of angles formula: l1*s1 + l2*(s1*c2 + c1*s2) + l3*(s1*c23 + c1*s23) = y
            # Rearrange terms:                                             q_s*s1 + q_c*c1 = y
            q_s = l1 + l2*c2 + l3*c23
            q_c = l2*s2 + l3*s23
            # Rearrange terms:                                                 q_s*s1 = y - q_c*c1
            # Square both sides:                                         q_s**2*s1**2 = y**2 - 2*y*q_c*c1 + q_c**2*c1**2
            # Rearrange terms    (q_c**2 + q_s**2)*c1**2 - 2*y*q_c*c1 + y**2 - q_s**2 = 0
            # Solve quadratic equation for c1
            a = q_c**2 + q_s**2
            b = -2*y*q_c
            c = y**2 - q_s**2

            discriminant = b**2 - 4*a*c
            if abs(discriminant) < epsilon:
                discriminant = 0.0
            elif discriminant < -epsilon:
                continue
            c1_1 = (-b + np.sqrt(discriminant))/(2*a)
            c1_2 = (-b - np.sqrt(discriminant))/(2*a)
            theta1_1 = np.arccos(c1_1)
            theta1_2 = np.arccos(c1_2)

            for theta1 in [theta1_1, theta1_2]:
                is_success, force, torques = findMaxForce(theta1, theta2, theta3)
                if is_success and force[0] < -MIN_NORMAL_REACTION:
                    print("Found thetas: " + str((theta1, theta2, theta3)))
                    valid_thetas.append((theta1, theta2, theta3))
    print("valid_thetas = " + str(valid_thetas))
    print("# valid_thetas = " + str(len(valid_thetas)))

if __name__ == "__main__":
    solve(0.3)
