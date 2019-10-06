import pydrake
from pydrake.solvers.mathematicalprogram import MathematicalProgram, Solve
from pydrake.symbolic import Variable
import numpy as np
import math
from enum import Enum
import time
import pdb
from constants import *

import eom

use_symbolic_derivs = False
STATE_SIZE = 8
TORQUE_SIZE = 3

HUB_MOTOR_MAX_TORQUE = 20
MIN_NORMAL_REACTION = 4
DISCRIMINANT_EPSILON = 1e-10
DYNAMICS_EPSILON = 1e-3

I_w = 1.0/2.0*m_w*(w_r**2)

NUM_TIME_STEPS = 50
TIME_INTERVAL = 1.0/NUM_TIME_STEPS

def calc_theta1_dd(state, u, is_symbolic):
    return eom.calc_theta1_dd(
            state[0], state[1], state[2],
            state[4], state[5], state[6],
            u[0], u[1], u[2], u[3]*w_r,
            is_symbolic)

def calc_theta2_dd(state, u, is_symbolic):
    return eom.calc_theta2_dd(
            state[0], state[1], state[2],
            state[4], state[5], state[6],
            u[0], u[1], u[2], u[3]*w_r,
            is_symbolic)

def calc_theta3_dd(state, u, is_symbolic):
    return eom.calc_theta3_dd(
            state[0], state[1], state[2],
            state[4], state[5], state[6],
            u[0], u[1], u[2], u[3]*w_r,
            is_symbolic)

def derivs(state, tau234, is_symbolic = True):
    if is_symbolic:
        sin = pydrake.symbolic.sin
        cos = pydrake.symbolic.cos
    else:
        sin = np.sin
        cos = np.cos

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

    tau1 = 0.0 # Unactuated
    tau = np.array([tau1, tau2, tau3, tau4])
    state_d = np.zeros_like(state)
    state_d[0:4] = state[4:8]
    state_d[4] = calc_theta1_dd(state, tau, is_symbolic)
    state_d[5] = calc_theta2_dd(state, tau, is_symbolic)
    state_d[6] = calc_theta3_dd(state, tau, is_symbolic)
    state_d[7] = tau234[2] / I_w
    return state_d

if __name__ == "__main__":
    mp = MathematicalProgram()
    state_over_time = np.zeros(shape=(NUM_TIME_STEPS, STATE_SIZE), dtype=pydrake.symbolic.Variable)
    tau234_over_time = np.zeros(shape=(NUM_TIME_STEPS, TORQUE_SIZE), dtype=pydrake.symbolic.Variable)

    initial_state = mp.NewContinuousVariables(8, "state_0")

    # Constrain initial velocity to be 0
    for j in range(4, 8):
        mp.AddConstraint(initial_state[j] <= 0.0)
        mp.AddConstraint(initial_state[j] >= 0.0)

    # Constrain initial theta to be between 0 ~ 180
    for j in range(0, 4):
        mp.AddConstraint(initial_state[j] >= 0.0)
        mp.AddConstraint(initial_state[j] <= np.pi)

    # for j in range(0, 3):
        # mp.AddConstraint(initial_state[j] <= np.pi/2.0)
        # mp.AddConstraint(initial_state[j] >= np.pi/2.0)

    # Constrain initial theta4 (i.e. front wheel) to be 0.0
    # mp.AddConstraint(initial_state[3] <= 0.0)
    # mp.AddConstraint(initial_state[3] >= 0.0)

    # Constrain initial front wheel position y position to be 0.0
    initial_wheel_position = eom.findFrontWheelPosition(initial_state[0], initial_state[1], initial_state[2])
    mp.AddConstraint(initial_wheel_position[1] <= 0.0)
    mp.AddConstraint(initial_wheel_position[1] >= 0.0)

    state_over_time[0] = initial_state

    for i in range(NUM_TIME_STEPS-1):
        print("Adding constraints for t = " + str(i))

        tau234 = mp.NewContinuousVariables(3, "tau234_%d" % i)
        next_state = mp.NewContinuousVariables(8, "state_%d" % (i+1))

        theta1 = next_state[0]
        theta2 = next_state[1]
        theta3 = next_state[2]
        theta4 = next_state[3]
        theta1_d = next_state[4]
        theta2_d = next_state[5]
        theta3_d = next_state[6]
        theta4_d = next_state[7]
        tau1 = 0.0 # Unactuated
        tau2 = tau234[0]
        tau3 = tau234[1]
        tau4 = tau234[2]

        # Add end force constraint
        # eom.calc_end_force_from_torques(
                # theta1, theta2, theta3,
                # theta1_d, theta2_d, theta3_d,
                # tau1, tau2, tau3)
        # mp.AddConstraint(COEFF_FRICTION*J.dot(tau123)[0] <= -w_r*tau234[2]) # FIXME

        # Constrain no x motion of front wheel
        # mp.AddConstraint(eom.findFrontWheelPosition(next_state[0], next_state[1], next_state[2])[0] <= initial_wheel_position[0])
        # mp.AddConstraint(eom.findFrontWheelPosition(next_state[0], next_state[1], next_state[2])[0] >= initial_wheel_position[0])

        # for j in range(3): # Constrain theta1, theta2, theta3
            # mp.AddConstraint(next_state[j] >= 0.0)
            # mp.AddConstraint(next_state[j] <= np.pi)

        if use_symbolic_derivs:
            for j in range(8):
                mp.AddConstraint(next_state[j] <= (state_over_time[i] + TIME_INTERVAL*derivs(state_over_time[i], tau234))[j])
                mp.AddConstraint(next_state[j] >= (state_over_time[i] + TIME_INTERVAL*derivs(state_over_time[i], tau234))[j])
        else:
            def next_state_constraint(stacked):
                next_state = stacked[0:STATE_SIZE]
                current_state = stacked[STATE_SIZE:STATE_SIZE*2]
                current_tau234 = stacked[STATE_SIZE*2:STATE_SIZE*2+TORQUE_SIZE]
                desired_state = current_state + TIME_INTERVAL*derivs(current_state, current_tau234, is_symbolic = False)
                diff = desired_state - next_state
                ret = np.zeros_like(stacked)
                ret[0:STATE_SIZE] = diff
                return ret

            stacked = np.concatenate([next_state, state_over_time[i], tau234])
            bounds = np.ones(stacked.shape)*DYNAMICS_EPSILON
            mp.AddConstraint(next_state_constraint, -bounds, bounds, stacked)

        tau234_over_time[i] = tau234
        state_over_time[i+1] = next_state

    mp.AddCost(0.01 * tau234_over_time[:,0].dot(tau234_over_time[:,0]))
    mp.AddCost(0.01 * tau234_over_time[:,1].dot(tau234_over_time[:,1]))
    mp.AddCost(0.01 * tau234_over_time[:,2].dot(tau234_over_time[:,2]))
    # This cost is incorrect, different case for if front wheel is left / right of back wheel
    # mp.AddCost(-(state_over_time[:,0].dot(state_over_time[:,0])))

    final_state = state_over_time[-1]
    # Constrain final velocity to be 0
    for j in range(4, 8):
        mp.AddConstraint(final_state[j] <= 0.0)
        mp.AddConstraint(final_state[j] >= 0.0)

    # Constrain final front wheel position
    final_front_wheel_pos = eom.findFrontWheelPosition(final_state[0], final_state[1], final_state[2])
    mp.AddConstraint(final_front_wheel_pos[1] <= STEP_HEIGHT)
    mp.AddConstraint(final_front_wheel_pos[1] >= STEP_HEIGHT)

    print("Begin solving...")
    t = time.time()
    result = Solve(mp)
    print("Done solving in " + str(time.time() - t) + "s")
    is_success = result.is_success()
    print("is_success = " + str(is_success))
    torque_over_time_star = result.GetSolution(tau234_over_time)
    state_over_time_star = result.GetSolution(state_over_time)

    for state in state_over_time_star:
        visualize.visualize(state[0], state[1], state[2])
        pdb.set_trace()
