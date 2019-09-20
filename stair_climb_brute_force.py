from pydrake.solvers.mathematicalprogram import MathematicalProgram, Solve
import numpy as np
from numpy import sin, cos
import math

HUB_MOTOR_STALL_TORQUE = 20.0
HUB_MOTOR_CONTINUOUS_STALL_TORQUE = HUB_MOTOR_STALL_TORQUE / 3.0
LINK_MOTOR_STALL_TORQUE = 5.0
LINK_MOTOR_CONTINUOUS_STALL_TORQUE = LINK_MOTOR_STALL_TORQUE / 3.0
COEFFICIENT_OF_FRICTION = 0.6
epsilon = 1e-10
g = 9.81
l_1 = 0.8
l_2 = 0.32
l_3 = l_1
l_b = 0.6*l_1 # length of battery link
m_2 = 0.2 # mass of l_1 - l_2 link motor
m_m = 0.4 # mass of battery link motor
m_b = 4.0 # mass of battery
m_3 = m_2 # mass of l_2 - l_3 motor
m_w = 4.0 # mass of wheel
w_r = 0.1 # wheel radius
HUB_MOTOR_FORCE = HUB_MOTOR_CONTINUOUS_STALL_TORQUE*w_r

def findMaxForce(theta1, theta2, theta3):
    theta12 = theta2 + theta1
    theta123 = theta3 + theta12
    s1 = sin(theta1)
    s2 = sin(theta2)
    s3 = sin(theta3)
    s12 = sin(theta12)
    s123 = sin(theta123)
    c1 = cos(theta1)
    c2 = cos(theta2)
    c3 = cos(theta3)
    c12 = cos(theta12)
    c123 = cos(theta123)

    J = np.array([
        [-l_1*s1 - l_2*s12 - l_3*s123, -l_2*s12 - l_3*s123, -l_3*s123],
        [l_1*c1 + l_2*c12 + l_3*c123, l_2*c12 + l_3*c123, l_3*c123]])

    tau1 = (-l_1*c1*m_2*g # torque due to l_1 - l_2 link motor
            + (-l_1*c1 + l_2/2.0*c12)*m_m*g # torque due to battery link motor
            + (-l_1*c1 + l_2/2.0*c12 + l_b*cos(theta12 + np.pi/2.0))*m_b*g # torque due to battery
            + (-l_1*c1 + l_2)*m_3*g # torque due to l_2 - l_3 link motor
            + (-l_1*c1 + l_2 + l_3*c3)*m_w*g # torque due to front wheel
    )

    # print("tau1 = " + str(tau1))

    mp = MathematicalProgram()

    tau23 = mp.NewContinuousVariables(2, "tau23")
    tau123 = np.append(tau1, tau23)
    mp.AddCost(J.dot(tau123)[0])
    mp.AddConstraint(tau23[0]**2 <= LINK_MOTOR_CONTINUOUS_STALL_TORQUE**2)
    mp.AddConstraint(tau23[1]**2 <= LINK_MOTOR_CONTINUOUS_STALL_TORQUE**2)
    mp.AddLinearConstraint(J.dot(tau123)[1] >= -HUB_MOTOR_FORCE/2.0)

    result = Solve(mp)
    is_success = result.is_success()
    torques = result.GetSolution(tau23)
    full_tau = np.append(tau1, torques)
    output_force = J.dot(full_tau)
    # print("is_success = " + str(is_success))
    # print("tau = " + str(full_tau))
    # print("F = " + str(J.dot(full_tau)))
    return is_success, output_force, torques

def solve(y, starting_theta2 = None, starting_theta3 = None):
    max_theta1 = 0.0
    best_theta123 = None
    best_torques = None
    if starting_theta2:
        theta2_range = np.arange(starting_theta2 - 0.1, starting_theta2 + 0.1, 0.02) # roughly 1 degree intervals
    else:
        theta2_range = np.arange(0.0, np.pi/2.0, 0.02) # roughly 1 degree intervals
    if starting_theta3:
        theta3_range = np.arange(starting_theta3 - 0.1, starting_theta3 + 0.1, 0.02)
    else:
        theta3_range = np.arange(0.0, np.pi/2.0, 0.02) # roughly 1 degree intervals
    for theta2 in theta2_range:
        for theta3 in theta3_range:
            s2 = sin(theta2)
            s3 = sin(theta3)
            s23 = sin(theta2+theta3)
            c2 = cos(theta2)
            c3 = cos(theta3)
            c23 = cos(theta2+theta3)
            # Given theta2, theta3 and y, find theta1
            # i.e. solve                                          l_1*s1 + l_2*s12 + l_3*s123 = y
            # Use sum of angles formula: l_1*s1 + l_2*(s1*c2 + c1*s2) + l_3*(s1*c23 + c1*s23) = y
            # Rearrange terms:                                                q_s*s1 + q_c*c1 = y
            q_s = l_1 + l_2*c2 + l_3*c23
            q_c = l_2*s2 + l_3*s23
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
                if is_success and force[0] < -(force[1]+HUB_MOTOR_FORCE)*COEFFICIENT_OF_FRICTION:
                    # print("Found thetas: " + str((theta1, theta2, theta3)))
                    # valid_thetas.append((theta1, theta2, theta3))
                    if theta1 > max_theta1:
                        max_theta1 = theta1
                        best_theta123 = (theta1, theta2, theta3)
                        best_torques = torques
                        print("new best theta123 = " + str(best_theta123))

    if best_theta123:
        print("thetas = " + str(best_theta123))
        print("torques = " + str(best_torques))
    else:
        print("No solution!")

if __name__ == "__main__":
    solve(0.1)
    solve(0.2)
    solve(0.3)
    solve(0.4)
    solve(0.5)
