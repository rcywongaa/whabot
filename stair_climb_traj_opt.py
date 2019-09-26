import pydrake
from pydrake.solvers.mathematicalprogram import MathematicalProgram, Solve
from pydrake.symbolic import (sin, cos, Variable)
import numpy as np
# from numpy import sin, cos
import math
from enum import Enum
import time
import pdb

import eom

STAIR_HEIGHT = 0.3
HUB_MOTOR_MAX_TORQUE = 20
MIN_NORMAL_REACTION = 4
EPSILON = 1e-10
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
I_w = 1.0/2.0*m_w*(w_r**2)
NUM_TIME_STEPS = 20
TIME_INTERVAL = 0.01
COEFF_FRICTION = 0.6

def calc_theta1_dd(state, u):
    return eom.calc_theta1_dd(
            state[0], state[1], state[2],
            state[4], state[5], state[6],
            u[0], u[1], u[2], u[3]*w_r)

def calc_theta2_dd(state, u):
    return eom.calc_theta2_dd(
            state[0], state[1], state[2],
            state[4], state[5], state[6],
            u[0], u[1], u[2], u[3]*w_r)

def calc_theta3_dd(state, u):
    return eom.calc_theta3_dd(
            state[0], state[1], state[2],
            state[4], state[5], state[6],
            u[0], u[1], u[2], u[3]*w_r)

def findTheta1(theta2, theta3, theta4):
    y = theta4*w_r
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
    if abs(discriminant) < EPSILON:
        discriminant = 0.0
    # elif discriminant < -EPSILON:
        # continue
    c1_1 = (-b + np.sqrt(discriminant))/(2*a)
    c1_2 = (-b - np.sqrt(discriminant))/(2*a)
    theta1_1 = np.arccos(c1_1)
    theta1_2 = np.arccos(c1_2)
    return theta1_1, theta1_2

def findTau1(theta1, theta2, theta3):
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

    tau1 = (-l_1*c1*m_2*g # torque due to l_1 - l_2 link motor
            + (-l_1*c1 + l_2/2.0*c12)*m_m*g # torque due to battery link motor
            + (-l_1*c1 + l_2/2.0*c12 + l_b*cos(theta12 + np.pi/2.0))*m_b*g # torque due to battery
            + (-l_1*c1 + l_2)*m_3*g # torque due to l_2 - l_3 link motor
            + (-l_1*c1 + l_2 + l_3*c3)*m_w*g # torque due to front wheel
    )
    return tau1

def find_I_wrt_1(theta1, theta2, theta3):
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

    I_2 = m_2*(l_1**2)
    I_3 = m_3*((l_1*s1 + l_2*s12)**2 + (l_1*c1 + l_2*c12)**2)
    I_4 = m_w*((l_1*s1 + l_2*s12 + l_3*s123)**2 + (l_1*c1 + l_2*c12 + l_3*c123)**2)
    return I_2 + I_3 + I_4

def find_I_wrt_2(theta2, theta3):
    s2 = sin(theta2)
    s3 = sin(theta3)
    s23 = sin(theta2 + theta3)
    c2 = cos(theta2)
    c3 = cos(theta3)
    c23 = cos(theta2 + theta3)

    I_3 = m_3*(l_2**2)
    I_4 = m_w*((l_2*s2 + l_3*s23)**2 + (l_2*c2 + l_3*c23)**2)
    return I_3 + I_4

def find_I_wrt_3():
    return m_w*(l_3**2)

def derivs(state, tau234):
    tau2 = tau234[0]
    tau3 = tau234[1]
    tau4 = tau234[2]
    theta1 = state[0]
    theta2 = state[1]
    theta3 = state[2]
    theta4 = state[3]
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

    tau1 = findTau1(theta1, theta2, theta3)
    tau = np.array([tau1, tau2, tau3, tau4])
    state_d = np.zeros_like(state)
    state_d[0:4] = state[4:8]
    state_d[4] = calc_theta1_dd(state, tau)
    state_d[5] = calc_theta2_dd(state, tau)
    state_d[6] = calc_theta3_dd(state, tau)
    state_d[7] = tau234[2] / I_w
    return state_d

def findJacobian(theta1, theta2, theta3):
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
    return J

if __name__ == "__main__":
    x = 0.92
    mp = MathematicalProgram()
    state_over_time = np.zeros(shape=(NUM_TIME_STEPS, 8), dtype=pydrake.symbolic.Expression)

    state_over_time[0] = mp.NewContinuousVariables(8, "state_0")

    tau234_over_time = np.zeros(shape=(NUM_TIME_STEPS, 3), dtype=pydrake.symbolic.Variable)
    for i in range(NUM_TIME_STEPS-1):
        print("Adding constraints for t = " + str(i))
        tau234 = mp.NewContinuousVariables(3, "tau234_%d" % i)
        tau234_over_time[i] = tau234

        state_over_time[i+1] = mp.NewContinuousVariables(8, "state_%d" % (i+1))

        theta1 = state_over_time[i+1][0]
        theta2 = state_over_time[i+1][1]
        theta3 = state_over_time[i+1][2]
        tau1 = findTau1(theta1, theta2, theta3)
        J = findJacobian(theta1, theta2, theta3)
        tau123 = np.array([tau1, tau234[0], tau234[2]])
        mp.AddConstraint(COEFF_FRICTION*J.dot(tau123)[0] <= -w_r*tau234[2])

        for j in range(3): # Constrain theta1, theta2, theta3
            mp.AddConstraint(state_over_time[i+1][j] >= 0.0)
            mp.AddConstraint(state_over_time[i+1][j] <= np.pi)
        for j in range(8):
            mp.AddConstraint(state_over_time[i+1][j] <= (state_over_time[i] + TIME_INTERVAL*derivs(state_over_time[i], tau234_over_time[i]))[j])
            mp.AddConstraint(state_over_time[i+1][j] >= (state_over_time[i] + TIME_INTERVAL*derivs(state_over_time[i], tau234_over_time[i]))[j])

    mp.AddCost(0.01 * tau234_over_time[:,0].dot(tau234_over_time[:,0]))
    mp.AddCost(0.01 * tau234_over_time[:,1].dot(tau234_over_time[:,1]))
    mp.AddCost(0.01 * tau234_over_time[:,2].dot(tau234_over_time[:,2]))
    mp.AddCost(-(state_over_time[:,0].dot(state_over_time[:,0])))
    target_theta4 = STAIR_HEIGHT / w_r

    # Constraint initial and final velocity to be 0
    for j in range(4, 8):
        mp.AddConstraint(state_over_time[0, j] <= 0.0)
        mp.AddConstraint(state_over_time[0, j] >= 0.0)
        mp.AddConstraint(state_over_time[-1, j] <= 0.0)
        mp.AddConstraint(state_over_time[-1, j] >= 0.0)

    for j in range(0, 4):
        mp.AddConstraint(state_over_time[0, j] >= 0.0)
        mp.AddConstraint(state_over_time[0, j] <= np.pi)

    mp.AddConstraint(state_over_time[0, 3] <= 0.0)
    mp.AddConstraint(state_over_time[0, 3] >= 0.0)

    mp.AddConstraint(state_over_time[-1, 3] <= target_theta4)
    mp.AddConstraint(state_over_time[-1, 3] >= target_theta4)

    print("Begin solving...")
    t = time.time()
    result = Solve(mp)
    print("Done solving in " + str(time.time() - t) + "s")
    is_success = result.is_success()
    print("is_success = " + str(is_success))
    torque_over_time = result.GetSolution(tau234_over_time)
    state_over_time = result.GetSolution(state_over_time)
    pdb.set_trace()
